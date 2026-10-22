!
! A superposition-of-atomic-densities guess from a JSON file, so a terco
! run can start where an mqc run would without mqc in the process.
!
! The file holds one density block per element for one basis, in that
! element's function order (shells as trc_basis_json builds them), Cartesian
! components always LEXICOGRAPHIC: mqc has no notion of TRC_CART_ORDER, so
! the format carries no field for it. trc_cart_perm permutes each block into
! this build's order on the way in. The guess for a molecule is those blocks
! on the diagonal in atom order. One file per basis:
!
!   {
!     "basis": "cc-pvdz",
!     "elements": {
!       "1": {"nao": 5,  "density": [ ... 25 numbers, column-major ... ]},
!       "8": {"nao": 15, "density": [ ... 225 numbers ... ]}
!     }
!   }
!
! The density is the TOTAL (alpha + beta) density; an unrestricted caller
! halves it per spin. A block whose size does not match the basis's count for
! that element is refused: it means a different basis and a silently wrong start.
!
! Normalisation: mqc's S_ii is 1 for s and p, 4 pi/5 for d and 4 pi/7 for f;
! terco's raw S_ii is 1 for all (`common_fac_sp` is cancelled by the basis
! builders, and is 1 for l >= 2). So s and p need no conversion, and a block
! element is divided by f(l_i) f(l_j), f(l) = sqrt((2l+1)/(4 pi)) /
! common_fac_sp(l), which is 1 for s and p.
!
! History: this module once divided by sqrt((2l+1)/(4 pi)) alone, which was only
! right because trc_basis_json double-folded `common_fac_sp` (two matching
! bugs); the divide was added when that was fixed. Energies cannot see a
! function's scale but a density can: Tr(D S) was 1.43 electrons for water and
! 10.0046 for cc-pVDZ with every energy check passing, which is why scf_sad
! integrates the guess against S.
!
module trc_guess_json
   use pic_types, only: dp
   use trc_constants, only: PI
   use trc_api, only: trc_basis_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_cart_perm, only: build_ao_perm
   use json_module, only: json_file
   implicit none
   private

   public :: trc_guess_from_json

contains

   !
   ! `dguess` (nao, nao) comes back allocated: the block-diagonal total density
   ! of the atoms of `b`. `at_z` are the atomic numbers in basis order.
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

      ! Functions per atom, where each atom's block starts, and each function's scale.
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

      ! Lexicographic -> build order, global AO index to global AO index (see the
      ! module header). `first(ia) + i - 1` is the global lexicographic index for any
      ! TRC_CART_ORDER, since only within-shell order changes.
      allocate (perm_l2b(b%nao), perm_b2l(b%nao))
      call build_ao_perm(b%nshell, b%sh_l, b%sh_ao, b%nao, perm_l2b, perm_b2l)

      allocate (dguess(b%nao, b%nao))
      dguess = 0.0_dp
      do ia = 1, b%natm
         write (buf, "(i0)") at_z(ia)
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
               ! (i, j) are local, lexicographic; map to global build-order AO indices.
               gi = perm_l2b(first(ia) + i - 1)
               gj = perm_l2b(first(ia) + j - 1)
               dguess(gi, gj) = blk(k)/(fac(gi)*fac(gj))
            end do
         end do
      end do
      call json%destroy()
      if (error%has_error() .and. allocated(dguess)) deallocate (dguess)
   end subroutine trc_guess_from_json

   integer function atom_of(b, ish) result(res)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: ish
      integer :: ia
      res = 0
      do ia = 1, b%natm
         if (all(abs(b%sh_r(:, ish) - b%at_r(:, ia)) < 1.0e-10_dp)) then
            res = ia
            return
         end if
      end do
   end function atom_of

   pure function itoa(i) result(t)
      integer, intent(in) :: i
      character(len=:), allocatable :: t
      character(len=16) :: buf
      write (buf, "(i0)") i
      t = trim(buf)
   end function itoa

   pure real(dp) function common_fac_sp(l) result(res)
      !! Duplicated from trc_api.F90 (and trc_ao_norm.F90, trc_basis_json.F90,
      !! trc_batch.F90) rather than imported -- see those modules' own comments
      !! on the same duplication. Used only to cancel itself out of `fac` for s
      !! and p, per this module's own header.
      integer, intent(in) :: l
      select case (l)
      case (0)
      res = 0.282094791773878143_dp
      case (1)
      res = 0.488602511902919921_dp
      case default
      res = 1.0_dp
      end select
   end function common_fac_sp

end module trc_guess_json
