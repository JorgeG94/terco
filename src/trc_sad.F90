!! Guesses built from free-atom calculations: SAD, and SADQ for ions
module trc_sad
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
   !! One restricted, spin-averaged Hartree-Fock run per element in that
   !! element's own shells, with fractional occupations so the atom is
   !! spherical, and the atomic densities placed on the diagonal of the
   !! molecular one in the molecule's function order. Cached per element,
   !! so cholesterol runs four atomic SCFs, not seventy-four.
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
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private
   public :: trc_sad_build, trc_sadq_build, trc_atom_ranges

contains

   subroutine trc_sad_build(b, dguess, error, verbose)
      !! `dguess(nao, nao)`: the total density, block diagonal over atoms.
      type(trc_basis_t), intent(in) :: b
      real(dp), allocatable, intent(out) :: dguess(:, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose
      call build_atomic(b, dguess, error, verbose)
   end subroutine trc_sad_build

   subroutine trc_sadq_build(b, nelec, dguess, error, verbose, qatom)
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
      !> Per-atom charge, positive for electrons removed. Overrides the
      !> proportional split; must have one entry per atom.
      real(dp), intent(in), optional :: qatom(:)

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
      call build_atomic(b, dguess, error, verbose, nel=nel)
   end subroutine trc_sadq_build

   !
   ! SAD and SADQ differ only in what each free atom is asked to hold, so they
   ! are one routine. `nel` absent is SAD: the neutral atom, unrestricted,
   ! Hund's rule, spread over the degenerate frontier. `nel` present is SADQ:
   ! that many electrons on that atom, restricted with fractional
   ! occupations, which is what a non-integer count needs.
   !
   subroutine build_atomic(b, dguess, error, verbose, nel)
      type(trc_basis_t), intent(in) :: b
      real(dp), allocatable, intent(out) :: dguess(:, :)
      type(error_t), intent(inout) :: error
      logical,  intent(in), optional :: verbose
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
      logical :: talk, sadq, hund

      talk = .false.
      if (present(verbose)) talk = verbose
      sadq = present(nel)
      allocate (dguess(b%nao, b%nao))
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
         a0 = b%sh_ao(first(ia))
         a1 = b%sh_ao(last(ia)) + ncart(b%sh_l(last(ia))) - 1
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
            call trc_scf_run(ab, na, nb, opts, res)
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

   pure subroutine hund_split(z, na, nb)
      !! Alpha and beta counts of the free atom's ground state: closed shells
      !! paired, the open subshell filled by Hund's rule, one spin first.
      !! Subshell order and capacities are the aufbau ones; the spread over
      !! the degenerate frontier is what keeps the atom spherical.
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
