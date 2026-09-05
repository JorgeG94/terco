!
! A superposition-of-atomic-densities guess from a JSON file, so a terco
! run can start where an mqc run would start without mqc in the process.
!
! The file holds one density block per element for one basis set, in that
! element's own function order -- the shells as trc_basis_json builds them
! from the same basis file, Cartesian components in terco's order -- and
! the guess for a molecule is those blocks placed on the diagonal, atom by
! atom, in the molecule's atom order. mqc writes it from its atomic guess
! (spherically averaged, restricted), one file per basis:
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
! differently, and neither is uniform in l on its own. terco folds
! libcint's common factor into the coefficients, so its S_ii is 1/(4 pi)
! for s, 3/(4 pi) for p, and 1 from d up. mqc's S_ii is 1 for s and p,
! 4 pi/5 for d, 4 pi/7 for f (measured, both of them). The RATIO is what
! matters here and it is 4 pi/(2l+1) for every l, so a block element is
! divided by f(l_i) f(l_j) with f(l) = sqrt((2l+1)/(4 pi)) on the way in.
! Energies cannot see a function's scale; a density can: Tr(D S) came out
! at 1.43 electrons for water before the s and p factors were understood,
! and at 10.0046 for cc-pVDZ before the d one was, with every energy check
! passing both times. That is why scf_sad integrates the guess against S.
!
module trc_guess_json
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t
   use trc_error, only: error_t, ERROR_VALIDATION
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
      integer :: ia, ish, n, k, i, j
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
         fac(b%sh_ao(ish):b%sh_ao(ish) + n - 1) = sqrt(real(2*b%sh_l(ish) + 1, dp)/(4.0_dp*PI))
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
               dguess(first(ia) + i - 1, first(ia) + j - 1) = &
                  blk(k)/(fac(first(ia) + i - 1)*fac(first(ia) + j - 1))
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

end module trc_guess_json
