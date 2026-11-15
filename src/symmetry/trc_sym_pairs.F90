module trc_sym_pairs
   !! Shell-pair orbit tables of a point group: the host-side input of the
   !! skeleton Fock build.
   !!
   !! Every operation g of the group sends shell s of the working basis to
   !! shell `shmap(s, g)` and a shell pair {a, b} to the pair {ga, gb}. Pairs
   !! are unordered, so the image of a pair stored as (a, b) may well be
   !! stored as (b', a'); the lookup is symmetric. A pair that maps to itself
   !! (possibly with a sign on its AO components, possibly swapped) is a pair
   !! with a nontrivial stabilizer: it simply has fewer distinct images.
   !!
   !! What is built
   !!  - `shmap`: the shell map of the *working* basis (the decontracted one
   !!    the four-centre build uses). It is read off the AO map of a
   !!    `trc_sym_ao_t` built on that basis, so the shell-to-shell
   !!    correspondence is the one the AO signs belong to, and it is checked
   !!    against the atom permutation `geom%perm`.
   !!  - `img(p, g)`: the pair that operation g sends pair p to.
   !!  - `rep(p)`: the representative of the orbit of p, the image with the
   !!    smallest canonical key (max(i,j)*(max(i,j)-1)/2 + min(i,j)). The key
   !!    does not depend on the order of the pair list, so every rank, and
   !!    every reordering of the list, picks the same representative.
   !!  - `osz(p)`: the orbit size, the number of distinct images of p.
   !!  - `rep_list`: the representatives, in ascending canonical key.
   !!  - `q(p)`: the Schwarz bound of p, taken as the maximum over its orbit.
   !!
   !! Screening that is invariant by construction
   !!  - Q: the symmetry-equivalent Q values differ in the last bits, and the
   !!    significance cut `Q_ab Q_max > thresh` can then keep p and drop gp.
   !!    `sym_pairs_symmetrize_q` replaces the triangular Q array by its orbit
   !!    maximum *before* the pairs are selected and binned, so the set of
   !!    significant pairs is closed and the size bucket is shared by an
   !!    orbit. `sym_pairs_build` refuses a pair list that is not closed.
   !!  - dsh(I,J), the largest |D| over a shell block: `sym_pairs_symmetrize_dsh`
   !!    takes the orbit maximum (O(nshell^2 * nop)); do it after every density
   !!    update, before the density screen.
   !!
   !! Weights. The weight of a pair as a bra is its orbit size `osz`, and the
   !! sum of `osz` over the representatives of a set closed under the group
   !! is the size of the set. Two further factors belong to the *plan*, not
   !! to this module, because they depend on which bin pair is launched:
   !!  - 1/2 on a diagonal bin pair (A, A), launched as the rectangle
   !!    reps(A) x A instead of the triangle;
   !!  - 1/2 on the quartets with ab = cd, which the existing digestion
   !!    already applies per quartet.
   !!
   !! How S14.2 consumes this (bins, plan, launch scalars)
   !!  1. Before `bins_build_binned_pairs`: call `sym_pairs_shell_map` once,
   !!     `sym_pairs_symmetrize_q` on the triangular Q. Then bin as today.
   !!  2. Call `sym_pairs_build` on the binned pair list (`sp_i`, `sp_j`,
   !!     `sp_q`; the arrays may be in any order).
   !!  3. Close the bins: the bin key of pair p is the key of `rep(p)` (the
   !!     size bucket, angular momenta and contraction degree are equal over an
   !!     orbit; the check tests the bucket). Within a bin store the
   !!     representatives first, ordered by `osz`, so each weight is a
   !!     contiguous sub-range. Split large bins only at orbit boundaries (the
   !!     cap becomes BIN_CAP + |G|).
   !!  4. For an unordered bin pair {A, B}, A /= B, launch the rectangle
   !!     reps(A) x B (or reps(B) x A, whichever has fewer quartets), once per
   !!     weight `osz`, with jfac and kfac multiplied by that weight. For the
   !!     diagonal pair (A, A) launch reps(A) x A with an extra 1/2.
   !!  5. Symmetrize the resulting skeleton G once, after the fold, with
   !!     `trc_sym_scf`'s symmetrizer.
   !!  Per bin pair the weighted count of quartets then equals the full one;
   !!  `check_sym_pairs` tests exactly that without any integral.
   !!
   !! Host code only. Nothing here is called by the SCF or the Fock build yet.
   use pic_types, only: dp, int64, int_index
   use pic_sorting, only: sort_index
   use pic_strings, only: to_string
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_basis_types, only: trc_basis_t
   use trc_cart, only: cart_ncart_of
   use trc_sym_ao, only: trc_sym_ao_t
   use trc_sym_geom, only: trc_sym_geom_t
   implicit none
   private

   public :: trc_sym_pairs_t, sym_pairs_build, sym_pairs_shell_map
   public :: sym_pairs_symmetrize_q, sym_pairs_symmetrize_dsh, sym_pairs_tri, sym_pairs_wanted

   type :: trc_sym_pairs_t
      !! Orbit tables of a pair list under a point group.
      integer :: nop = 0  !! number of operations
      integer :: nshell = 0  !! shells of the working basis
      integer :: npair = 0  !! pairs in the list
      integer :: nrep = 0  !! number of orbits
      integer, allocatable :: shmap(:, :)  !! (nshell, nop) image shell
      integer, allocatable :: img(:, :)  !! (npair, nop) image pair
      integer, allocatable :: rep(:)  !! (npair) representative of the orbit
      integer, allocatable :: osz(:)  !! (npair) orbit size
      integer, allocatable :: rep_list(:)  !! (nrep) representatives, ascending canonical key
      real(dp), allocatable :: q(:)  !! (npair) Schwarz bound, maximum over the orbit
   contains
      procedure :: size_counts => pairs_size_counts
   end type trc_sym_pairs_t

contains

   pure function sym_pairs_tri(i, j) result(k)
      !! Canonical triangular key of the unordered pair {i, j}.
      integer, intent(in) :: i, j
      integer(int64) :: k
      integer(int64) :: hi, lo
      hi = max(i, j)
      lo = min(i, j)
      k = hi*(hi - 1)/2 + lo
   end function sym_pairs_tri

   subroutine sym_pairs_shell_map(bas, ao, geom, shmap, error)
      !! shmap(s, g): the shell that operation g sends working-basis shell s
      !! to. `ao` must have been built on `bas` (`sym_ao_build`). Refuses a
      !! map that is not a permutation of the shells, that does not carry
      !! the first AO of a shell to the first AO of a shell, or that moves a
      !! shell to an atom other than `geom%perm` says.
      type(trc_basis_t), intent(in) :: bas
      type(trc_sym_ao_t), intent(in) :: ao
      type(trc_sym_geom_t), intent(in) :: geom
      integer, allocatable, intent(out) :: shmap(:, :)
      type(error_t), intent(inout) :: error

      integer, allocatable :: first(:), shell_of(:)
      logical, allocatable :: hit(:)
      integer :: s, g, nsh, nao, t, c

      nsh = bas%nshell
      nao = ao%nao
      if (.not. allocated(ao%ao_map) .or. ao%nop < 1 .or. nsh < 1) then
         call error%set(ERROR_VALIDATION, "sym_pairs: the AO maps are not built")
         return
      end if
      if (.not. allocated(geom%perm) .or. size(geom%perm, 2) /= ao%nop .or. size(geom%perm, 1) /= bas%natm) then
         call error%set(ERROR_VALIDATION, "sym_pairs: the atom permutation does not match the basis and the AO maps")
         return
      end if
      if (nao /= merge(bas%nao_sph, bas%nao, ao%spherical) .or. (ao%spherical .neqv. bas%spherical)) then
         call error%set(ERROR_VALIDATION, "sym_pairs: the AO maps were not built on this basis")
         return
      end if

      allocate (first(nsh), shell_of(nao), shmap(nsh, ao%nop), hit(nsh))
      shell_of = 0
      do s = 1, nsh
         first(s) = merge(bas%sh_ao_sph(s), bas%sh_ao(s), ao%spherical)
         c = merge(2*bas%sh_l(s) + 1, cart_ncart_of(bas%sh_l(s)), ao%spherical)
         if (first(s) < 1 .or. first(s) + c - 1 > nao) then
            call error%set(ERROR_VALIDATION, "sym_pairs: shell "//to_string(s)//" lies outside the AO space")
            return
         end if
         shell_of(first(s):first(s) + c - 1) = s
      end do

      do g = 1, ao%nop
         hit = .false.
         do s = 1, nsh
            t = shell_of(ao%ao_map(first(s), g))
            if (t < 1 .or. ao%ao_map(first(s), g) /= first(t)) then
               call error%set(ERROR_VALIDATION, "sym_pairs: operation "//to_string(g)//" does not map shell "// &
                              to_string(s)//" onto a whole shell")
               return
            end if
            if (bas%sh_at(t) /= geom%perm(bas%sh_at(s), g) .or. bas%sh_l(t) /= bas%sh_l(s)) then
               call error%set(ERROR_VALIDATION, "sym_pairs: the shell map of operation "//to_string(g)// &
                              " disagrees with the atom permutation")
               return
            end if
            if (hit(t)) then
               call error%set(ERROR_VALIDATION, "sym_pairs: operation "//to_string(g)//" is not a permutation of the shells")
               return
            end if
            hit(t) = .true.
            shmap(s, g) = t
         end do
      end do
   end subroutine sym_pairs_shell_map

   subroutine sym_pairs_symmetrize_q(shmap, q)
      !! Replace the triangular Schwarz array q(sym_pairs_tri(i, j)) by its
      !! maximum over the orbit of {i, j}. Idempotent. Call before the
      !! significant pairs are selected and binned.
      integer, intent(in) :: shmap(:, :)
      real(dp), intent(inout) :: q(:)
      real(dp), allocatable :: qn(:)
      integer :: i, j, g
      integer(int64) :: k

      if (size(shmap, 2) <= 1) return
      allocate (qn(size(q)))
      qn = q
      do i = 1, size(shmap, 1)
         do j = 1, i
            k = sym_pairs_tri(i, j)
            do g = 2, size(shmap, 2)
               qn(k) = max(qn(k), q(sym_pairs_tri(shmap(i, g), shmap(j, g))))
            end do
         end do
      end do
      q = qn
   end subroutine sym_pairs_symmetrize_q

   subroutine sym_pairs_symmetrize_dsh(shmap, dsh)
      !! Replace the shell-block density maxima dsh(I, J) by the maximum over
      !! the orbit of (I, J). Idempotent.
      integer, intent(in) :: shmap(:, :)
      real(dp), intent(inout) :: dsh(:, :)
      real(dp), allocatable :: dn(:, :)
      integer :: i, j, g

      if (size(shmap, 2) <= 1) return
      allocate (dn(size(dsh, 1), size(dsh, 2)))
      dn = dsh
      do j = 1, size(shmap, 1)
         do i = 1, size(shmap, 1)
            do g = 2, size(shmap, 2)
               dn(i, j) = max(dn(i, j), dsh(shmap(i, g), shmap(j, g)))
            end do
         end do
      end do
      dsh = dn
   end subroutine sym_pairs_symmetrize_dsh

   function sym_pairs_wanted(shmap, sh_l) result(want)
      !! want(i, j): the ordered shell pair (i, j) is one the unique-pair build
      !! must form primitive pairs for. Those are the orbit representatives
      !! (the member of an orbit with the smallest canonical key) in the
      !! orders that are used:
      !!  - (i, j) with i >= j, the order the Schwarz bound is computed in;
      !!  - the order of every image: the bins store a pair with its higher-l
      !!    shell first (index order breaks a tie), and the image of a pair
      !!    (a, b) stored so is (ga, gb), in this order, whichever of the
      !!    orbit's members it is. An image is rebuilt from the representative
      !!    in that order by sign flips alone, so that order of the
      !!    representative has to exist.
      !! Needs no pair list, only the shell map and the shells' l: O(nshell^2
      !! nop) integer work. The ordered pairs of all other orbit members are
      !! rebuilt, never formed.
      integer, intent(in) :: shmap(:, :), sh_l(:)
      logical, allocatable :: want(:, :)
      logical, allocatable :: isrep(:, :)
      integer :: i, j, g, nsh, nop, a, b, u, v
      integer(int64) :: k0, k

      nsh = size(shmap, 1)
      nop = size(shmap, 2)
      allocate (want(nsh, nsh), isrep(nsh, nsh))
      want = .false.
      do j = 1, nsh
         do i = 1, nsh
            k0 = sym_pairs_tri(i, j)
            isrep(i, j) = .true.
            do g = 2, nop
               k = sym_pairs_tri(shmap(i, g), shmap(j, g))
               if (k < k0) then
                  isrep(i, j) = .false.
                  exit
               end if
            end do
         end do
      end do
      do i = 1, nsh
         do j = 1, i
            if (isrep(i, j)) want(i, j) = .true.
            if (sh_l(j) > sh_l(i)) then
               a = j
               b = i
            else
               a = i
               b = j
            end if
            do g = 1, nop
               u = shmap(a, g)
               v = shmap(b, g)
               if (isrep(u, v)) want(u, v) = .true.
            end do
         end do
      end do
   end function sym_pairs_wanted

   function sym_pairs_build(shmap, sp_i, sp_j, sp_q, error) result(tab)
      !! Orbit tables of the pair list (sp_i(p), sp_j(p)) with Schwarz bounds
      !! sp_q(p), under the shell map `shmap`. The list must hold each
      !! unordered pair at most once and be closed under the group; otherwise
      !! `error` is set (symmetrize Q before selecting the pairs).
      integer, intent(in) :: shmap(:, :)
      integer, intent(in) :: sp_i(:), sp_j(:)
      real(dp), intent(in) :: sp_q(:)
      type(error_t), intent(inout) :: error
      type(trc_sym_pairs_t) :: tab

      integer, allocatable :: pidx(:)
      integer(int64) :: key, kb
      integer :: nsh, nop, np, p, g, m, i, j, best, cnt

      nsh = size(shmap, 1)
      nop = size(shmap, 2)
      np = size(sp_i)
      if (nop < 1 .or. nsh < 1 .or. size(sp_j) /= np .or. size(sp_q) /= np) then
         call error%set(ERROR_VALIDATION, "sym_pairs_build: inconsistent argument sizes")
         return
      end if
      if (nsh > 60000) then
         call error%set(ERROR_VALIDATION, "sym_pairs_build: more shells than the pair keys can number")
         return
      end if

      allocate (pidx(int(sym_pairs_tri(nsh, nsh))))
      pidx = 0
      do p = 1, np
         i = sp_i(p)
         j = sp_j(p)
         if (min(i, j) < 1 .or. max(i, j) > nsh) then
            call error%set(ERROR_VALIDATION, "sym_pairs_build: pair "//to_string(p)//" names a shell that does not exist")
            return
         end if
         key = sym_pairs_tri(i, j)
         if (pidx(key) /= 0) then
            call error%set(ERROR_VALIDATION, "sym_pairs_build: the pair ("//to_string(i)//","//to_string(j)// &
                           ") is listed twice")
            return
         end if
         pidx(key) = p
      end do

      tab%nop = nop
      tab%nshell = nsh
      tab%npair = np
      tab%shmap = shmap
      allocate (tab%img(np, nop), tab%rep(np), tab%osz(np), tab%q(np))
      do p = 1, np
         do g = 1, nop
            m = pidx(sym_pairs_tri(shmap(sp_i(p), g), shmap(sp_j(p), g)))
            if (m == 0) then
               call error%set(ERROR_VALIDATION, "sym_pairs_build: the image of pair ("//to_string(sp_i(p))//","// &
                              to_string(sp_j(p))//") under operation "//to_string(g)// &
                              " is not in the list; the pair set is not closed (symmetrize Q first)")
               return
            end if
            tab%img(p, g) = m
         end do
      end do

      do p = 1, np
         best = p
         kb = sym_pairs_tri(sp_i(p), sp_j(p))
         cnt = 0
         tab%q(p) = 0.0_dp
         do g = 1, nop
            m = tab%img(p, g)
            tab%q(p) = max(tab%q(p), sp_q(m))
            key = sym_pairs_tri(sp_i(m), sp_j(m))
            if (key < kb) then
               kb = key
               best = m
            end if
            ! distinct images only; an orbit has at most nop members
            if (.not. any(tab%img(p, 1:g - 1) == m)) cnt = cnt + 1
         end do
         tab%rep(p) = best
         tab%osz(p) = cnt
      end do

      tab%nrep = count(tab%rep == [(p, p=1, np)])
      allocate (tab%rep_list(tab%nrep))
      m = 0
      do p = 1, np
         if (tab%rep(p) == p) then
            m = m + 1
            tab%rep_list(m) = p
         end if
      end do
      ! ascending canonical key
      block
         integer(int64), allocatable :: kk(:)
         integer(int_index), allocatable :: ord(:)
         integer :: a
         allocate (kk(tab%nrep), ord(tab%nrep))
         do a = 1, tab%nrep
            kk(a) = sym_pairs_tri(sp_i(tab%rep_list(a)), sp_j(tab%rep_list(a)))
         end do
         ! ascending pair key; the keys of distinct pairs are distinct
         call sort_index(kk, ord)
         tab%rep_list = tab%rep_list(ord)
      end block
   end function sym_pairs_build

   function pairs_size_counts(this) result(cnt)
      !! cnt(n): the number of orbits of size n (n = 1..nop).
      class(trc_sym_pairs_t), intent(in) :: this
      integer :: cnt(max(this%nop, 1))
      integer :: k, p
      cnt = 0
      do k = 1, this%nrep
         p = this%rep_list(k)
         cnt(this%osz(p)) = cnt(this%osz(p)) + 1
      end do
   end function pairs_size_counts

end module trc_sym_pairs
