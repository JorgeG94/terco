!
! Runtime selection of Rys over HGP, per class, from TRC_RYS_CLASSES.
!
! Up to HGP_LMAX both schemes have a kernel, and the default is per class:
! the classes in gen_rys.py's FOCK_RYS_DEFAULT (rys_default_class) go to the
! chunked Rys kernels, where that measured faster than HGP on every system
! tried, and the rest to HGP. TRC_RYS_CLASSES overrides the default for the
! whole build. Above HGP_LMAX there is no HGP kernel at all -- HGP stops at
! d -- so a class with any shell above it goes to Rys whatever
! TRC_RYS_CLASSES says. This module is the only place that decides which one
! a launch uses, and it is hand-written rather than generated because the
! decision itself -- parse an environment variable, look up a key -- has
! nothing to do with angular momentum and nothing to regenerate if LMAX or
! TRC_CART_ORDER changes; the default list itself is generated with the
! kernels it names.
!
! READ ONCE PER FOCK BUILD, ON THE HOST. `rys_select_init` is an ordinary
! subroutine called from trc_eri.F90 before any kernel launches; the result
! is a plain logical table that `rys_select_route` turns into one branch at each
! launch site in trc_binkernel.F90, before a `do concurrent` starts. Nothing
! here runs on the device, and no kernel's own do-concurrent-only rule is
! touched by it.
!
module trc_rys_select
   use pic_types, only: dp
   use trc_rys_kernels, only: rys_has_class, rys_default_class
   use trc_rys_table1, only: rys_table_init1
   use trc_rys_table2, only: rys_table_init2
   use trc_rys_table3, only: rys_table_init3
   use trc_rys_table4, only: rys_table_init4
   use trc_rys_table5, only: rys_table_init5
   use trc_rys_table6, only: rys_table_init6
   use trc_rys_table7, only: rys_table_init7
   use trc_rys_table8, only: rys_table_init8
   use trc_rys_table9, only: rys_table_init9
   implicit none
   private
   public :: rys_select_init, rys_select_route, rys_select_tables_init

   integer, parameter :: RYS_CLASS_RADIX = 11
   !! gen_perclass.py's CLASS_RADIX, fixed at 11 so a dispatch key means the
   !! same thing at every LMAX -- see its own comment there. Duplicated as a
   !! literal here rather than imported, because importing it would pull in
   !! trc_pc_kernels (HGP's dispatcher module) for one integer this module
   !! only needs to parse an environment variable. The two must agree; if
   !! CLASS_RADIX ever changes, update this one too.
   integer, parameter :: RYS_KEY_MAX = RYS_CLASS_RADIX**4 - 1

   integer, parameter :: HGP_LMAX = 2
   !! Highest shell angular momentum with an HGP per-class kernel. The HGP
   !! set is generated at --lmax 2 and kept at that ceiling at every
   !! TERCO_LMAX: its d kernels generated at --lmax 3 are identical, and its
   !! f kernels would be 175 classes of HGP text for integrals Rys builds
   !! anyway.

   logical, save :: rys_all = .false.
   logical, save :: rys_default = .true.
   !! TRC_RYS_CLASSES unset, empty or "default": route by rys_default_class.
   logical, save :: rys_table(0:RYS_KEY_MAX) = .false.

contains

   subroutine rys_select_tables_init()
      !! Push the Rys root tables (1..9 roots) to the device, as boys_init does
      !! the Boys table. Their `declare create` allocates device storage; the
      !! tables' own comment says it does not carry the host initialiser.
      !! nvfortran evidently initialises the device copy anyway (GPU Rys
      !! energies matched the host build to 1e-12 before anything called
      !! this), but that is undocumented compiler behaviour, not a contract.
      !! Called by every launcher of a kernel that reads the tables, before
      !! the launch: trc_eri's build (Fock, Schwarz and gradient kernels),
      !! trc_integrals' trc_df_2c/trc_df_3c and trc_scf_grad_df (the DF
      !! derivative kernels). Repeating it is harmless: one small `update
      !! device` per table, onto whichever device is current.
      call rys_table_init1()
      call rys_table_init2()
      call rys_table_init3()
      call rys_table_init4()
      call rys_table_init5()
      call rys_table_init6()
      call rys_table_init7()
      call rys_table_init8()
      call rys_table_init9()
   end subroutine rys_select_tables_init

   subroutine rys_select_init()
      !! Unset, empty or "default" keeps the measured per-class default
      !! (rys_default_class); "none" sends every class up to HGP_LMAX to HGP;
      !! "all" every class the build generated a kernel for to Rys; a
      !! comma-separated list of four-digit class tags (e.g. "2222,2221,1111")
      !! selects exactly those that also have a kernel, and HGP for the rest.
      !! A tag that is not exactly four ASCII digits is ignored rather than
      !! failing the whole value over one typo.
      character(len=256) :: value
      integer :: stat, length, p, q, key
      rys_all = .false.
      rys_default = .true.
      rys_table = .false.
      call get_environment_variable("TRC_RYS_CLASSES", value, length, stat)
      if (stat /= 0 .or. length == 0) return
      if (trim(value) == "default") return
      rys_default = .false.
      if (trim(value) == "all") then
         rys_all = .true.
         return
      end if
      if (trim(value) == "none") return
      p = 1
      do while (p <= length)
         q = index(value(p:length), ",")
         if (q == 0) then
            q = length
         else
            q = p + q - 2
         end if
         if (q >= p) then
            key = class_key(value(p:q))
            if (key >= 0) rys_table(key) = .true.
         end if
         p = q + 2
      end do
   end subroutine rys_select_init

   pure logical function rys_select_route(key) result(res)
      !! Should THIS launch -- class key `key`, the same fixed-radix key
      !! trc_binkernel.F90 forms from (sLA,sLB,sLC,sLD) -- go to Rys? Folds in
      !! `rys_has_class`, so a class TRC_RYS_CLASSES asks for but that has no
      !! generated kernel falls back to HGP instead of reaching an unbuilt
      !! reference or, with "all", silently skipping the class entirely. A
      !! class above HGP_LMAX always goes to Rys, and has no fallback: HGP's
      !! dispatcher would run nothing for it, so a missing Rys kernel there
      !! is a stop, not a silently absent block of integrals.
      integer, intent(in) :: key
      logical :: wanted
      if (key_lmax(key) > HGP_LMAX) then
         if (.not. rys_has_class(key)) error stop &
            "terco: a class above d has no Rys kernel in this build"
         res = .true.
         return
      end if
      if (rys_all) then
         wanted = .true.
      else if (rys_default) then
         wanted = rys_default_class(key)
      else if (key >= 0 .and. key <= RYS_KEY_MAX) then
         wanted = rys_table(key)
      else
         wanted = .false.
      end if
      res = wanted .and. rys_has_class(key)
   end function rys_select_route

   pure integer function key_lmax(key) result(res)
      !! The highest of the four shell angular momenta packed into `key`.
      integer, intent(in) :: key
      integer :: k, i
      k = key
      res = 0
      do i = 1, 4
         res = max(res, mod(k, RYS_CLASS_RADIX))
         k = k/RYS_CLASS_RADIX
      end do
   end function key_lmax

   pure integer function class_key(tag) result(res)
      !! Four ASCII digits "0".."9" -- la, lb, lc, ld in order -- to the same
      !! fixed-radix key gen_perclass.py's CLASS_RADIX forms. -1 if `tag` is
      !! not exactly four digits.
      character(len=*), intent(in) :: tag
      integer :: i, d
      res = 0
      if (len_trim(tag) /= 4) then
         res = -1
         return
      end if
      do i = 1, 4
         d = index("0123456789", tag(i:i)) - 1
         if (d < 0) then
            res = -1
            return
         end if
         res = res*RYS_CLASS_RADIX + d
      end do
   end function class_key

end module trc_rys_select
