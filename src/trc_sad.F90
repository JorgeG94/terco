!! The superposition-of-atomic-densities guess, built here
module trc_sad
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
   !! does; the first SCF iteration absorbs the difference.
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private
   public :: trc_sad_build

contains

   subroutine trc_sad_build(b, dguess, error, verbose)
      !! `dguess(nao, nao)`: the total density, block diagonal over atoms.
      type(trc_basis_t), intent(in) :: b
      real(dp), allocatable, intent(out) :: dguess(:, :)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose

      integer, allocatable :: atom_of(:), first(:), last(:)
      integer, allocatable :: cache_z(:), cache_n(:)
      type dblock
         real(dp), allocatable :: d(:, :)
      end type dblock
      type(dblock), allocatable :: cache(:)
      type(trc_basis_t) :: ab
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      integer :: ia, ish, ncached, k, a0, a1, z, nz, nel
      logical :: talk

      talk = .false.
      if (present(verbose)) talk = verbose
      allocate (dguess(b%nao, b%nao))
      dguess = 0.0_dp

      ! shell -> atom by centre, and each atom's contiguous shell range
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
            call error%set(ERROR_VALIDATION, "trc_sad: shell "//itoa(ish)//" sits on no atom")
            return
         end if
         ia = atom_of(ish)
         if (first(ia) == 0) first(ia) = ish
         if (last(ia) /= 0 .and. last(ia) /= ish - 1) then
            call error%set(ERROR_VALIDATION, "trc_sad: the shells of atom "//itoa(ia)//" are not contiguous")
            return
         end if
         last(ia) = ish
      end do

      allocate (cache_z(b%natm), cache_n(b%natm), cache(b%natm))
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
         k = 0
         do ish = 1, ncached
            if (cache_z(ish) == z .and. cache_n(ish) == nz) k = ish
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
            nel = z
            opts%nelec_frac = real(nel, dp)
            call trc_scf_run(ab, nel/2, nel/2, opts, res)
            if (.not. res%converged) then
               call error%set(ERROR_VALIDATION, "trc_sad: the atomic SCF for Z = "//itoa(z)//" did not converge: "// &
                              trim(res%message))
               call ab%release()
               return
            end if
            if (talk) print '(a,i0,a,i0,a,f16.8,a,i0,a)', "  sad: Z = ", z, " (", nz, " functions) E = ", &
               res%energy, " in ", res%iterations, " iterations"
            ncached = ncached + 1
            cache_z(ncached) = z; cache_n(ncached) = nz
            allocate (cache(ncached)%d(nz, nz))
            cache(ncached)%d = res%dmat(:, :, 1)
            call ab%release()
            k = ncached
         end if
         dguess(a0:a1, a0:a1) = cache(k)%d
      end do
   end subroutine trc_sad_build

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
