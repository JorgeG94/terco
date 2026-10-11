!
! A superposition-of-atomic-densities guess from a JSON file, so a terco
! run can start where an mqc run would start without mqc in the process.
!
! The file holds one density block per element for one basis set, in that
! element's own function order -- the shells as trc_basis_json builds them
! from the same basis file, Cartesian components in LEXICOGRAPHIC order, a
! FIXED CONVENTION: mqc always writes lexicographic, regardless of which
! order this build was built with (mqc has no notion of TRC_CART_ORDER),
! and there is no field in the format to say so because it never varies.
! trc_cart_perm permutes each block's components into THIS build's own
! order on the way in, below, so every caller downstream sees an ordinary
! build-order density -- the guess for a molecule is those blocks placed on
! the diagonal, atom by atom, in the molecule's atom order. mqc writes it
! from its atomic guess (spherically averaged, restricted), one file per
! basis:
!
!   {
!     "basis": "cc-pvdz",
!     "elements": {
!       "1": {"nao": 5,  "density": [ ... 25 numbers, column-major ... ]},
!       "8": {"nao": 15, "density": [ ... 225 numbers ... ]}
!     }
!   }
!
! The density is the TOTAL (alpha + beta) density: it is what a restricted
! SCF starts from, and an unrestricted caller halves it per spin. A block
! whose size does not match the basis's count for that element is refused,
! since a mismatch means a different basis and a silently wrong start.
!
! CONVENTION. The two codes scale a Cartesian function x^l e^{-a r^2}
! differently for d and up, and neither is uniform in l on its own. mqc's
! S_ii is 1 for s and p, 4 pi/5 for d, 4 pi/7 for f (measured). terco's own
! raw S_ii is 1 for s and p too -- `trc_basis_t%build` folds libcint's
! `common_fac_sp` into an s or p shell's coefficients, but every builder
! that feeds it (this file's `trc_basis_json`, test/trc_test_basis.F90's
! `build_631g`) hands it a coefficient already divided by that factor, so
! the two cancel and terco's raw <s|s>/<p_x|p_x> come out 1, matching mqc
! (and pyscf/libcint) directly: no conversion needed for s or p. `build`
! does not touch d and up (common_fac_sp(l>=2) = 1), so terco's OWN d/f
! stays at 1 while mqc's is 4 pi/5, 4 pi/7 -- a REAL difference, needing a
! REAL conversion there. The ratio is 4 pi/(2l+1), so a block element is
! divided by f(l_i) f(l_j) with f(l) = sqrt((2l+1)/(4 pi)) / common_fac_sp(l)
! on the way in -- the common_fac_sp divide makes f(l) = 1 for s and p
! (nothing to convert) and leaves it sqrt((2l+1)/(4 pi)) for d and up
! (common_fac_sp(l>=2) = 1, a no-op divide there).
!
! THIS MODULE USED TO DIVIDE BY sqrt((2l+1)/(4 pi)) ALONE, s and p included
! -- which happened to be the right number only because trc_basis_json had
! its OWN bug of the same shape (folding common_fac_sp twice into an s or p
! shell, making terco's raw S_ii for those come out 1/(4 pi) and 3/(4 pi)
! for real, not 1). Two bugs the same shape, canceling: that is also why
! Tr(D S) came out right (10.0046 for cc-pVDZ, see below) despite both. The
! common_fac_sp divide was added here the moment trc_basis_json's fold was
! fixed and this module's own conversion stopped canceling it. Energies
! cannot see a function's scale; a density can: Tr(D S) came out at 1.43
! electrons for water before the s and p factors were understood, and at
! 10.0046 for cc-pVDZ before the d one was, with every energy check passing
! both times. That is why scf_sad integrates the guess against S.
!
module trc_guess_json
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_cart_perm, only: build_ao_perm
   use json_module, only: json_file
   implicit none
   private

   public :: trc_guess_from_json

   real(dp), parameter :: PI = 3.14159265358979323846_dp

contains

   !
   ! `dguess` (nao, nao) comes back allocated, the block-diagonal total
   ! density of the atoms of `b`. `at_z` are the atomic numbers in the order
   ! the basis was built; `b%sh_r` against `b%at_r` says which shells belong
   ! to which atom.
   !
   subroutine trc_guess_from_json(path, b, at_z, dguess, error)
      character(len=*), intent(in) :: path
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: at_z(b%natm)
      real(dp), allocatable, intent(out) :: dguess(:, :)
      type(error_t), intent(inout) :: error

      type(json_file) :: json
      character(len=:), allocatable :: key
      real(dp), allocatable :: blk(:), fac(:)
      integer, allocatable :: nfun(:), first(:)
      integer, allocatable :: perm_l2b(:), perm_b2l(:)
      integer :: ia, ish, n, k, i, j, gi, gj
      logical :: found
      character(len=16) :: buf

      ! Functions per atom and where each atom's block starts, from the
      ! shell centres: shells are built atom by atom in atom order. And the
      ! scale of every function, for the convention above.
      allocate (nfun(b%natm), first(b%natm), fac(b%nao))
      nfun = 0
      do ish = 1, b%nshell
         ia = atom_of(b, ish)
         n = (b%sh_l(ish) + 1)*(b%sh_l(ish) + 2)/2
         nfun(ia) = nfun(ia) + n
         fac(b%sh_ao(ish):b%sh_ao(ish) + n - 1) = sqrt(real(2*b%sh_l(ish) + 1, dp)/(4.0_dp*PI)) &
                                                 /common_fac_sp(b%sh_l(ish))
      end do
      first(1) = 1
      do ia = 2, b%natm
         first(ia) = first(ia - 1) + nfun(ia - 1)
      end do
      if (first(b%natm) + nfun(b%natm) - 1 /= b%nao) then
         call error%set(ERROR_VALIDATION, "trc_guess_json: shells are not grouped by atom in atom order")
         return
      end if

      call json%initialize()
      call json%load(filename=path)
      if (json%failed()) then
         call error%set(ERROR_VALIDATION, "trc_guess_json: cannot read "//trim(path))
         call json%destroy()
         return
      end if

      ! LEXICOGRAPHIC -> this build's order, global AO index -> global AO
      ! index (see the module header). `first(ia) + i - 1` for a LOCAL,
      ! lexicographic-order `i` (mqc's own convention, 1-based) is the
      ! correct GLOBAL LEXICOGRAPHIC AO index regardless of TRC_CART_ORDER:
      ! shell order by atom does not change, only the within-shell component
      ! order does, and `first`/`nfun` above are built from counts alone.
      allocate (perm_l2b(b%nao), perm_b2l(b%nao))
      call build_ao_perm(b%nshell, b%sh_l, b%sh_ao, b%nao, perm_l2b, perm_b2l)

      allocate (dguess(b%nao, b%nao))
      dguess = 0.0_dp
      do ia = 1, b%natm
         write (buf, '(i0)') at_z(ia)
         key = "elements."//trim(buf)
         call json%get(key//".nao", n, found)
         if (.not. found) then
            call error%set(ERROR_VALIDATION, "trc_guess_json: no block for element "//trim(buf)//" in "//trim(path))
            exit
         end if
         if (n /= nfun(ia)) then
            call error%set(ERROR_VALIDATION, "trc_guess_json: element "//trim(buf)//" has "// &
                           itoa(n)//" functions in "//trim(path)//" and "//itoa(nfun(ia))// &
                           " in this basis: a different basis set")
            exit
         end if
         call json%get(key//".density", blk, found)
         if (.not. found .or. size(blk) /= n*n) then
            call error%set(ERROR_VALIDATION, "trc_guess_json: malformed density block for element "//trim(buf))
            exit
         end if
         k = 0
         do j = 1, n
            do i = 1, n
               k = k + 1
               ! (i, j) are LOCAL, LEXICOGRAPHIC (mqc's file). Map to this
               ! build's global AO indices before placing the element, so
               ! `dguess` comes back ordinary build-order data like every
               ! other AO-indexed quantity this library produces.
               gi = perm_l2b(first(ia) + i - 1)
               gj = perm_l2b(first(ia) + j - 1)
               dguess(gi, gj) = blk(k)/(fac(gi)*fac(gj))
            end do
         end do
      end do
      call json%destroy()
      if (error%has_error() .and. allocated(dguess)) deallocate (dguess)
   end subroutine trc_guess_from_json

   integer function atom_of(b, ish)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: ish
      integer :: ia
      atom_of = 0
      do ia = 1, b%natm
         if (all(abs(b%sh_r(:, ish) - b%at_r(:, ia)) < 1.0e-10_dp)) then
            atom_of = ia
            return
         end if
      end do
   end function atom_of

   pure function itoa(i) result(t)
      integer, intent(in) :: i
      character(len=:), allocatable :: t
      character(len=16) :: buf
      write (buf, '(i0)') i
      t = trim(buf)
   end function itoa

   pure real(dp) function common_fac_sp(l)
      !! Duplicated from trc_api.F90 (and trc_ao_norm.F90, trc_basis_json.F90,
      !! trc_batch.F90) rather than imported -- see those modules' own comments
      !! on the same duplication. Used only to cancel itself out of `fac` for s
      !! and p, per this module's own header.
      integer, intent(in) :: l
      select case (l)
      case (0); common_fac_sp = 0.282094791773878143_dp
      case (1); common_fac_sp = 0.488602511902919921_dp
      case default; common_fac_sp = 1.0_dp
      end select
   end function common_fac_sp

end module trc_guess_json
