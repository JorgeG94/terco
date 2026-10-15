module trc_sad
   !! Guesses built from free-atom calculations: SAD, and SADQ for ions
   !!
   !! Two entry points over the same machinery.
   !!
   !! `trc_sad_build` is SAD, the superposition of atomic densities, and is
   !! described below.
   !!
   !! `trc_sadq_build` is SADQ, the same construction with the atoms CHARGED:
   !! the molecule's charge is spread over its atoms and each one's density
   !! comes from a free-atom SCF at its own fractional electron count. SAD
   !! hands a cation the neutral atoms' density and lets the first iteration
   !! absorb the difference, which for a highly charged system is a long way
   !! to come back from. On a neutral molecule the two agree by construction.
   !!
   !! SAD's own description, which SADQ inherits except for the occupations:
   !!
   !! One Hartree-Fock run per element in that element's own shells --
   !! UNRESTRICTED with integer, Hund's-rule occupations for a neutral
   !! atom asking for its own Z (see `build_atomic`'s "WHICH ATOM TO RUN"
   !! comment: spin-restricting oxygen costs 0.124 Ha and an extra SCF
   !! iteration for water), restricted with fractional occupations for a
   !! SADQ atom carrying a non-integer charge -- then SPHERICALLY AVERAGED
   !! (`average_atom_density`, below) before it is placed on the diagonal
   !! of the molecular density, in the molecule's function order. An open
   !! Hund shell (O's 2p, say) is not itself spherical: its minority-spin
   !! electron sits in whichever orientation the eigensolver's degenerate
   !! subspace happened to return, which is not reproducible across GPUs,
   !! compilers or cuSOLVER versions, only the atom's total ENERGY is.
   !! Averaging removes that orientation after the fact, keeping the lower
   !! Hund energy and the atom's electron count exactly, so the guess
   !! (and the SCF path it starts) is deterministic again. Cached per
   !! element, so cholesterol runs four atomic SCFs, not seventy-four.
   !!
   !! The atomic basis is `subset` of the molecular one over the atom's shell
   !! range, with the atom list cut down to that one nucleus; `subset` copies
   !! the folded coefficients as they are, which is the point of using it.
   !! The shells of one atom are assumed contiguous, which every reader here
   !! produces; a basis that interleaves atoms is refused.
   !!
   !! Charged molecules get the neutral atoms' density unscaled, as PySCF
   !! does; the first SCF iteration absorbs the difference. That sentence is
   !! why SADQ exists -- it is the one thing SAD will not do for an ion.
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_linalg, only: trc_linalg_t
   use trc_sph_tables, only: SPH_T
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private
   public :: trc_sad_build, trc_sadq_build, trc_atom_ranges, average_atom_density

contains

   subroutine trc_sad_build(b, dguess, error, verbose, la)
      !! `dguess(nao, nao)`: the total density, block diagonal over atoms.
      type(trc_basis_t), intent(in) :: b
      real(dp), allocatable, intent(out) :: dguess(:, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose
      type(trc_linalg_t), intent(inout), optional :: la
      !! the calculation's cuBLAS/cuSOLVER handles, for the atomic SCFs
      call build_atomic(b, dguess, error, verbose, la=la)
   end subroutine trc_sad_build

   subroutine trc_sadq_build(b, nelec, dguess, error, verbose, qatom, la)
      !! SADQ: SAD with the atoms carrying the molecule's charge.
      !!
      !! `nelec` is the MOLECULE's electron count, so the charge it has to
      !! place is sum(Z) - nelec. Without `qatom` that charge goes onto the
      !! atoms in proportion to Z, which is the only parameter-free split
      !! there is: it gets the total electron count and a monotone ordering
      !! right and claims nothing about chemistry. It is not a charge model
      !! -- a cation's electron does not leave uniformly, it leaves the
      !! least electronegative atom -- so `qatom` is there for a caller who
      !! knows better, and a caller who does should use it.
      !!
      !! On a neutral molecule this is SAD, up to the occupations: the
      !! atomic SCF runs restricted with fractional occupations rather than
      !! unrestricted with Hund's rule, because a fractional electron count
      !! has no integer per-spin split to give it.
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: nelec
      real(dp), allocatable, intent(out) :: dguess(:, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose
      type(trc_linalg_t), intent(inout), optional :: la
      !! the calculation's cuBLAS/cuSOLVER handles, for the atomic SCFs
      real(dp), intent(in), optional :: qatom(:)
      !! Per-atom charge, positive for electrons removed. Overrides the
      !! proportional split; must have one entry per atom.

      real(dp), allocatable :: nel(:)
      real(dp) :: ztot
      integer :: ia

      ztot = sum(b%at_z(1:b%natm))
      allocate (nel(b%natm))
      if (present(qatom)) then
         if (size(qatom) /= b%natm) then
            call error%set(ERROR_VALIDATION, "trc_sadq: qatom has "//trim(itoa(size(qatom)))// &
                           " entries for "//trim(itoa(b%natm))//" atoms")
            return
         end if
         nel = b%at_z(1:b%natm) - qatom
      else if (ztot > 0.0_dp) then
         nel = b%at_z(1:b%natm)*(real(nelec, dp)/ztot)
      else
         nel = 0.0_dp
      end if
      do ia = 1, b%natm
         if (nel(ia) < 0.0_dp) then
            call error%set(ERROR_VALIDATION, "trc_sadq: atom "//trim(itoa(ia))// &
                           " is left with a negative electron count")
            return
         end if
      end do
      call build_atomic(b, dguess, error, verbose, nel=nel, la=la)
   end subroutine trc_sadq_build

   !
   ! SAD and SADQ differ only in what each free atom is asked to hold, so they
   ! are one routine. `nel` absent is SAD: the neutral atom, unrestricted,
   ! Hund's rule, spread over the degenerate frontier. `nel` present is SADQ:
   ! that many electrons on that atom, restricted with fractional
   ! occupations, which is what a non-integer count needs.
   !
   subroutine build_atomic(b, dguess, error, verbose, nel, la)
      type(trc_basis_t), intent(in) :: b
      real(dp), allocatable, intent(out) :: dguess(:, :)
      type(error_t), intent(inout) :: error
      logical,  intent(in), optional :: verbose
      type(trc_linalg_t), intent(inout), optional :: la
      !! the calculation's cuBLAS/cuSOLVER handles, for the atomic SCFs
      real(dp), intent(in), optional :: nel(:)

      integer, allocatable :: atom_of(:), first(:), last(:)
      integer, allocatable :: cache_z(:), cache_n(:)
      real(dp), allocatable :: cache_e(:)
      type dblock
         real(dp), allocatable :: d(:, :)
      end type dblock
      type(dblock), allocatable :: cache(:)
      type(trc_basis_t) :: ab
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      integer :: ia, ish, ncached, k, a0, a1, z, nz, nel_i, na, nb
      real(dp) :: want
      logical :: talk, sadq, hund, sph

      talk = .false.
      if (present(verbose)) talk = verbose
      sadq = present(nel)
      ! `dguess` is in whatever AO space `b` itself exposes to the SCF
      ! (trc_api.F90's `spherical`/`nao_sph`): a spherical molecular basis
      ! gets a spherical SAD guess, built from spherical ATOMIC SCFs below,
      ! because `b%subset` (and so every `ab` built from it) carries the
      ! same `spherical` flag -- `trc_scf_run(ab, ...)` on a spherical atom
      ! already returns `res%dmat` at `ab%nao_sph`, not at the atom's
      ! CARTESIAN function count, so the block this routine places into
      ! `dguess` has to be addressed and sized to match.
      sph = b%spherical
      allocate (dguess(merge(b%nao_sph, b%nao, sph), merge(b%nao_sph, b%nao, sph)))
      dguess = 0.0_dp

      call trc_atom_ranges(b, atom_of, first, last, error)
      if (error%has_error()) return

      allocate (cache_z(b%natm), cache_n(b%natm), cache_e(b%natm), cache(b%natm))
      ncached = 0
      opts%functional = ""
      opts%guess = "core"
      opts%frac_occ = .true.
      opts%conv_energy = 1.0e-8_dp
      opts%conv_diis = 1.0e-5_dp
      opts%max_iter = 100
      opts%verbose = .false.

      do ia = 1, b%natm
         z = nint(b%at_z(ia))
         if (first(ia) == 0) cycle   ! a ghost, or an atom with no functions
         if (sph) then
            a0 = b%sh_ao_sph(first(ia))
            a1 = b%sh_ao_sph(last(ia)) + (2*b%sh_l(last(ia)) + 1) - 1
         else
            a0 = b%sh_ao(first(ia))
            a1 = b%sh_ao(last(ia)) + ncart(b%sh_l(last(ia))) - 1
         end if
         nz = a1 - a0 + 1
         want = real(z, dp)
         if (sadq) want = nel(ia)
         ! The electron count joins the cache key: two atoms of one element
         ! share an SCF only if they are asked for the same charge, which
         ! under the proportional split they always are.
         k = 0
         do ish = 1, ncached
            if (cache_z(ish) == z .and. cache_n(ish) == nz .and. &
                abs(cache_e(ish) - want) < 1.0e-12_dp) k = ish
         end do
         if (k == 0) then
            call b%subset(first(ia), last(ia), ab)
            ! the isolated atom: one nucleus, this one
            ab%natm = 1
            deallocate (ab%at_z, ab%at_r)
            allocate (ab%at_z(1), ab%at_r(3, 1))
            ab%at_z(1) = b%at_z(ia)
            ab%at_r(:, 1) = b%at_r(:, ia)
            call ab%to_device()
            opts%nelec_frac = want
            !
            ! WHICH ATOM TO RUN, DECIDED PER ATOM.
            !
            ! An atom asked for its own Z has an integer count and an
            ! integer per-spin split, so it gets the unrestricted Hund atom
            ! -- which is a BETTER atom: spin-restricting oxygen costs it
            ! 0.124 Ha, and starting water from the restricted atoms cost an
            ! extra SCF iteration. Only a fractional count needs the
            ! restricted fractional path, and it is the only one that has
            ! it. So SADQ on a neutral molecule is SAD exactly, rather than
            ! approximately, and differs only where it has something to say.
            !
            hund = abs(want - real(z, dp)) < 1.0e-12_dp
            opts%unrestricted = hund
            if (hund) then
               call hund_split(z, na, nb)
            else
               ! nspin is 1 here, and `nocc` is taken from nelec_frac, so
               ! these two only have to be legal.
               nel_i = int(want/2.0_dp)
               na = nel_i; nb = nel_i
            end if
            call trc_scf_run(ab, na, nb, opts, res, la=la)
            if (.not. res%converged) then
               call error%set(ERROR_VALIDATION, merge("trc_sadq", "trc_sad ", sadq)// &
                              ": the atomic SCF for Z = "//trim(itoa(z))//" did not converge: "// &
                              trim(res%message))
               call ab%release()
               return
            end if
            if (talk) then
               if (sadq .and. .not. hund) then
                  print '(a,i0,a,i0,a,f10.6,a,f16.8,a,i0,a)', "  sadq: Z = ", z, " (", nz, &
                     " functions) n = ", want, " E = ", res%energy, " in ", res%iterations, " iterations"
               else
                  print '(a,i0,a,i0,a,f16.8,a,i0,a)', "  sad: Z = ", z, " (", nz, " functions) E = ", &
                     res%energy, " in ", res%iterations, " iterations"
               end if
            end if
            ncached = ncached + 1
            cache_z(ncached) = z; cache_n(ncached) = nz; cache_e(ncached) = want
            allocate (cache(ncached)%d(nz, nz))
            if (hund) then
               cache(ncached)%d = res%dmat(:, :, 1) + res%dmat(:, :, 2)
            else
               cache(ncached)%d = res%dmat(:, :, 1)
            end if
            ! SPHERICALLY AVERAGE. Run unconditionally -- a closed-shell or
            ! SADQ-fractional atom is already (to SCF convergence) exactly
            ! what this produces, so averaging it is a no-op; a Hund atom
            ! is not, and this is the one place that gets fixed. See
            ! `average_atom_density`'s own header for the derivation.
            call average_atom_density(ab, sph, cache(ncached)%d)
            call ab%release()
            k = ncached
         end if
         dguess(a0:a1, a0:a1) = cache(k)%d
      end do
   end subroutine build_atomic

   !
   ! shell -> atom by centre, and each atom's contiguous shell range.
   !
   subroutine trc_atom_ranges(b, atom_of, first, last, error)
      type(trc_basis_t), intent(in) :: b
      integer, allocatable, intent(out) :: atom_of(:), first(:), last(:)
      type(error_t), intent(inout) :: error
      integer :: ish, ia

      allocate (atom_of(b%nshell), first(b%natm), last(b%natm))
      first = 0; last = 0
      do ish = 1, b%nshell
         atom_of(ish) = 0
         do ia = 1, b%natm
            if (all(abs(b%sh_r(:, ish) - b%at_r(:, ia)) < 1.0e-10_dp)) then
               atom_of(ish) = ia
               exit
            end if
         end do
         if (atom_of(ish) == 0) then
            call error%set(ERROR_VALIDATION, "trc_sad: shell "//trim(itoa(ish))//" sits on no atom")
            return
         end if
         ia = atom_of(ish)
         if (first(ia) == 0) first(ia) = ish
         if (last(ia) /= 0 .and. last(ia) /= ish - 1) then
            call error%set(ERROR_VALIDATION, "trc_sad: the shells of atom "//trim(itoa(ia))// &
                           " are not contiguous")
            return
         end if
         last(ia) = ish
      end do
   end subroutine trc_atom_ranges

   !
   ! SPHERICAL AVERAGING. `d`, the atom's own total (alpha+beta) density in
   ! `ab`'s AO space, is replaced in place by its rotational average over
   ! every shell pair: zero for an (l_i, l_j) pair with l_i /= l_j, and
   ! (trace/(2l+1)) times the identity for one with l_i == l_j.
   !
   ! WHY THIS IS THE RIGHT AVERAGE. Two shells of an isolated atom with the
   ! same l are two copies of the same SO(3) irreducible representation (the
   ! angular integral of two real solid harmonics factors out of their
   ! radial overlap and is a Kronecker delta in m, REGARDLESS of either
   ! shell's radial part), so the atom's own overlap S restricted to that
   ! shell pair is already exactly a scalar times the identity -- Schur's
   ! lemma, the only SO(3)-invariant map between two copies of one
   ! irreducible representation. Two shells with l_i /= l_j are
   ! INEQUIVALENT irreps and admit no invariant map at all: that block of S
   ! is exactly zero. A rotational average of ANY operator commutes with
   ! every rotation by construction, so the averaged density must look like
   ! S block by block: a multiple of the identity where l_i == l_j, zero
   ! elsewhere -- and the multiple that preserves the block's own trace
   ! (hence preserves Tr(D_block S_block) = (scalar)*Tr(D_block), hence the
   ! atom's electron count Tr(D S) exactly, summed over every block) is
   ! trace/(2l+1). A different-l block already contributes zero to Tr(D S)
   ! before this runs (S is zero there), so zeroing it costs nothing either.
   !
   ! CARTESIAN ATOMS, AND WHY THIS ADDS A DELTA RATHER THAN REPLACING `d`.
   ! A Cartesian shell's six d (or ten f, ...) components are not the
   ! orthonormal-per-m basis the identity above needs: six Cartesian d =
   ! five spherical d + an s-type CONTAMINANT (x^2+y^2+z^2), ten Cartesian
   ! f = seven spherical f + a p-type contaminant (r^2 x, r^2 y, r^2 z), and
   ! so on -- see trc_sph.F90's own header. T = SPH_T (the block trc_sph.F90
   ! applies to S/H/D/F) embeds only the PURE l component into the Cartesian
   ! shell, so
   ! T^T T is NOT the identity for l >= 2 (checked numerically against a
   ! real shell's own Cartesian overlap: T^T S T = I holds -- T IS correctly
   ! normalised against the atom's own metric -- but T^T T alone does not),
   ! and replacing `d` outright by T * average(T^T d T) * T^T would both
   ! DISCARD the contaminant's electron density (a real, non-zero
   ! contribution once a Cartesian f shell polarises into an occupied p
   ! orbital, say) and fail to reproduce `d` itself on the *clean* round
   ! trip T^T d T when `d` already has no contaminant -- either way, Tr(d S)
   ! would move.
   !
   ! The fix is to apply only the DELTA the averaging makes in the pure
   ! representation, through T, and leave the rest of `d` (whatever is not
   ! pure-spherical -- the contaminant, and anything T^T(.)T cannot see)
   ! untouched:
   !
   !     d := d + T * (average(T^T d T) - T^T d T) * T^T
   !
   ! This preserves Tr(d S) EXACTLY, for any T, any metric, any contaminant,
   ! because it only needs that S_sph := T^T S T is exactly zero for an
   ! (l_i, l_j) pair with l_i /= l_j and exactly s_ij*(identity) for one
   ! with l_i == l_j -- true by the same angular-orthogonality argument as
   ! above, applied to PURE embeddings only, regardless of either shell's
   ! contaminant -- and that the delta itself has EXACTLY ZERO trace on
   ! every (l_i, l_j) block with l_i == l_j (it is average(block) minus
   ! block, and average(block) was built to have the same trace) and is
   ! multiplied by a ZERO S_sph block whenever l_i /= l_j. Tr(delta * S_sph)
   ! is a sum of terms that are each either `scalar * 0` or `0 * Tr(delta
   ! block) = 0`, so it is exactly zero, and Tr((d + T delta T^T) S) =
   ! Tr(d S) + Tr(delta * S_sph) = Tr(d S). `cart_to_sph_atom`/
   ! `sph_to_cart_atom` are a plain host copy of trc_sph.F90's own T-block
   ! contraction (no `do concurrent`, no device data assumed -- this runs
   ! once per cached element on an atom-sized matrix, never inside a device
   ! region).
   !
   subroutine average_atom_density(ab, sph, d)
      type(trc_basis_t), intent(in) :: ab
      logical, intent(in) :: sph
      real(dp), intent(inout) :: d(:, :)
      real(dp), allocatable :: ds(:, :), ds_avg(:, :), ddelta(:, :)

      if (sph) then
         call average_sph_blocks(ab%nshell, ab%sh_l, ab%sh_ao_sph, ab%nao_sph, d)
      else
         allocate (ds(ab%nao_sph, ab%nao_sph), ds_avg(ab%nao_sph, ab%nao_sph))
         allocate (ddelta(ab%nao, ab%nao))
         call cart_to_sph_atom(ab, d, ds)
         ds_avg = ds
         call average_sph_blocks(ab%nshell, ab%sh_l, ab%sh_ao_sph, ab%nao_sph, ds_avg)
         call sph_to_cart_atom(ab, ds_avg - ds, ddelta)
         d = d + ddelta
         deallocate (ds, ds_avg, ddelta)
      end if
   end subroutine average_atom_density

   subroutine average_sph_blocks(nshell, sh_l, sh_a, n, d)
      !! The Schur's-lemma average itself, on a SPHERICAL (n, n) matrix:
      !! see `average_atom_density`'s header. `sh_a` is each shell's
      !! first-AO offset in `d` (1-based, m ascending) -- `ab%sh_ao_sph`
      !! from the caller, local to the atom.
      integer,  intent(in) :: nshell, n
      integer,  intent(in) :: sh_l(nshell), sh_a(nshell)
      real(dp), intent(inout) :: d(n, n)
      integer :: ish, jsh, li, lj, nsi, nsj, bi, bj, r
      real(dp) :: tr, avg

      do jsh = 1, nshell
         do ish = 1, nshell
            li = sh_l(ish); lj = sh_l(jsh)
            bi = sh_a(ish); bj = sh_a(jsh)
            nsi = 2*li + 1; nsj = 2*lj + 1
            if (li /= lj) then
               d(bi:bi + nsi - 1, bj:bj + nsj - 1) = 0.0_dp
            else
               tr = 0.0_dp
               do r = 0, nsi - 1
                  tr = tr + d(bi + r, bj + r)
               end do
               avg = tr/real(nsi, dp)
               d(bi:bi + nsi - 1, bj:bj + nsi - 1) = 0.0_dp
               do r = 0, nsi - 1
                  d(bi + r, bj + r) = avg
               end do
            end if
         end do
      end do
   end subroutine average_sph_blocks

   subroutine cart_to_sph_atom(ab, dc, ds)
      !! ds := T^T dc T, shell pair by shell pair -- the same contraction
      !! trc_sph.F90's `cart_to_sph_driver` runs, as a plain sequential
      !! host loop (no `do concurrent`, no device data assumed): see
      !! `average_atom_density`'s header for why this is its own copy
      !! rather than a call into that module.
      type(trc_basis_t), intent(in) :: ab
      real(dp), intent(in)  :: dc(ab%nao, ab%nao)
      real(dp), intent(out) :: ds(ab%nao_sph, ab%nao_sph)
      integer :: ish, jsh, li, lj, nci, ncj, nsi, nsj, ai, aj, bi, bj, p, q, r, s
      real(dp) :: acc, inner

      ds = 0.0_dp
      do jsh = 1, ab%nshell
         do ish = 1, ab%nshell
            li = ab%sh_l(ish); lj = ab%sh_l(jsh)
            nci = ncart(li); ncj = ncart(lj)
            nsi = 2*li + 1; nsj = 2*lj + 1
            ai = ab%sh_ao(ish); aj = ab%sh_ao(jsh)
            bi = ab%sh_ao_sph(ish); bj = ab%sh_ao_sph(jsh)
            do s = 1, nsj
               do r = 1, nsi
                  acc = 0.0_dp
                  do q = 1, ncj
                     inner = 0.0_dp
                     do p = 1, nci
                        inner = inner + SPH_T(p, r, li)*dc(ai + p - 1, aj + q - 1)
                     end do
                     acc = acc + inner*SPH_T(q, s, lj)
                  end do
                  ds(bi + r - 1, bj + s - 1) = acc
               end do
            end do
         end do
      end do
   end subroutine cart_to_sph_atom

   subroutine sph_to_cart_atom(ab, ds, dc)
      !! dc := T ds T^T, the reverse of `cart_to_sph_atom`; see its header.
      type(trc_basis_t), intent(in) :: ab
      real(dp), intent(in)  :: ds(ab%nao_sph, ab%nao_sph)
      real(dp), intent(out) :: dc(ab%nao, ab%nao)
      integer :: ish, jsh, li, lj, nci, ncj, nsi, nsj, ai, aj, bi, bj, p, q, r, s
      real(dp) :: acc, inner

      dc = 0.0_dp
      do jsh = 1, ab%nshell
         do ish = 1, ab%nshell
            li = ab%sh_l(ish); lj = ab%sh_l(jsh)
            nci = ncart(li); ncj = ncart(lj)
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
   end subroutine sph_to_cart_atom

   pure subroutine hund_split(z, na, nb)
      !! Alpha and beta counts of the free atom's ground state: closed shells
      !! paired, the open subshell filled by Hund's rule, one spin first.
      !! Subshell order and capacities are the aufbau ones -- the spread over
      !! the degenerate frontier gives the right trace per block but not yet
      !! a spherical atom (an open shell's degenerate occupation can land in
      !! any orientation); `average_atom_density`, above, is what makes the
      !! atom (and so the guess) spherical and orientation-independent.
      integer, intent(in) :: z
      integer, intent(out) :: na, nb
      integer, parameter :: cap(19) = [2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 6]
      integer :: left, i, half, part
      left = z; na = 0; nb = 0
      do i = 1, size(cap)
         if (left <= 0) exit
         part = min(left, cap(i))
         half = cap(i)/2
         if (part <= half) then
            na = na + part
         else
            na = na + half
            nb = nb + (part - half)
         end if
         left = left - part
      end do
      ! whatever is beyond the table, paired
      na = na + (left + 1)/2; nb = nb + left/2
   end subroutine hund_split

   pure integer function ncart(l)
      integer, intent(in) :: l
      ncart = (l + 1)*(l + 2)/2
   end function ncart

   pure function itoa(i) result(t)
      integer, intent(in) :: i
      character(len=12) :: t
      write (t, '(i0)') i
   end function itoa

end module trc_sad
