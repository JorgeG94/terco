!
! Runtime selection of Rys over HGP, per class, from TRC_RYS_CLASSES.
!
! Rys is a measurement tool next to HGP, not a replacement for it: this
! module is the only place that decides which one a launch uses, and it is
! hand-written rather than generated because the decision itself -- parse an
! environment variable, look up a key -- has nothing to do with angular
! momentum and nothing to regenerate if LMAX or TRC_CART_ORDER changes.
!
! READ ONCE PER FOCK BUILD, ON THE HOST. `rys_select_init` is an ordinary
! subroutine called from trc_fock.F90 before any kernel launches; the result
! is a plain logical table that `rys_route` turns into one branch at each
! launch site in trc_binkernel.F90, before a `do concurrent` starts. Nothing
! here runs on the device, and no kernel's own do-concurrent-only rule is
! touched by it.
!
module trc_rys_select
   use trc_boys, only: dp
   use trc_rys_kernels, only: rys_has_class
   implicit none
   private
   public :: rys_select_init, rys_route

   integer, parameter :: RYS_CLASS_RADIX = 11
   !! gen_perclass.py's CLASS_RADIX, fixed at 11 so a dispatch key means the
   !! same thing at every LMAX -- see its own comment there. Duplicated as a
   !! literal here rather than imported, because importing it would pull in
   !! trc_pc_kernels (HGP's dispatcher module) for one integer this module
   !! only needs to parse an environment variable. The two must agree; if
   !! CLASS_RADIX ever changes, update this one too.
   integer, parameter :: RYS_KEY_MAX = RYS_CLASS_RADIX**4 - 1

   logical, save :: rys_all = .false.
   logical, save :: rys_table(0:RYS_KEY_MAX) = .false.

contains

   subroutine rys_select_init()
      !! "none" (unset, or anything unrecognised) selects no class; "all"
      !! every class TERCO_RYS generated a kernel for; a comma-separated list
      !! of four-digit class tags (e.g. "2222,2221,1111") selects exactly
      !! those that also have a kernel. A tag that is not exactly four ASCII
      !! digits is ignored rather than failing the whole value over one typo.
      character(len=256) :: val
      integer :: stat, length, p, q, key
      rys_all = .false.
      rys_table = .false.
      call get_environment_variable("TRC_RYS_CLASSES", val, length, stat)
      if (stat /= 0 .or. length == 0) return
      if (trim(val) == "all") then
         rys_all = .true.
         return
      end if
      if (trim(val) == "none") return
      p = 1
      do while (p <= length)
         q = index(val(p:length), ",")
         if (q == 0) then
            q = length
         else
            q = p + q - 2
         end if
         if (q >= p) then
            key = class_key(val(p:q))
            if (key >= 0) rys_table(key) = .true.
         end if
         p = q + 2
      end do
   end subroutine rys_select_init

   pure logical function rys_route(key)
      !! Should THIS launch -- class key `key`, the same fixed-radix key
      !! trc_binkernel.F90 forms from (sLA,sLB,sLC,sLD) -- go to Rys? Folds in
      !! `rys_has_class`, so a class TRC_RYS_CLASSES asks for but that has no
      !! generated kernel falls back to HGP instead of reaching an unbuilt
      !! reference or, with "all", silently skipping the class entirely.
      integer, intent(in) :: key
      logical :: wanted
      if (rys_all) then
         wanted = .true.
      else if (key >= 0 .and. key <= RYS_KEY_MAX) then
         wanted = rys_table(key)
      else
         wanted = .false.
      end if
      rys_route = wanted .and. rys_has_class(key)
   end function rys_route

   pure integer function class_key(tag)
      !! Four ASCII digits "0".."9" -- la, lb, lc, ld in order -- to the same
      !! fixed-radix key gen_perclass.py's CLASS_RADIX forms. -1 if `tag` is
      !! not exactly four digits.
      character(len=*), intent(in) :: tag
      integer :: i, d
      class_key = 0
      if (len_trim(tag) /= 4) then
         class_key = -1
         return
      end if
      do i = 1, 4
         d = index("0123456789", tag(i:i)) - 1
         if (d < 0) then
            class_key = -1
            return
         end if
         class_key = class_key*RYS_CLASS_RADIX + d
      end do
   end function class_key

end module trc_rys_select
