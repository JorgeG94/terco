!
! A basis from a MolSSI Basis Set Exchange JSON file, the files metalquicha
! ships under basis_sets/, so a terco run can take the orbital and auxiliary
! sets mqc would hand it without mqc in the process.
!
! Ported from metalquicha's mqc_json_basis_reader, with the same reading:
! one shell per COEFFICIENT COLUMN, so an SP shell becomes an s and a p
! shell and a general contraction (cc-pVDZ oxygen: nine s primitives, three
! columns) becomes three shells over the same primitives -- which is the
! segmentation terco's kernels want and what the mqc bridge does before
! handing a basis over. Values are read one at a time by indexed path: BSE
! stores them as strings to keep every digit.
!
! Coefficients are normalised as metalquicha and pyscf normalise them: each
! primitive to unit self-overlap, then the contraction to unit self-overlap.
! The function_type ("gto_spherical" on a d shell) is read and ignored:
! terco is Cartesian, and a spherical set is run in its Cartesian form, as
! pyscf does with cart=True.
!
module trc_basis_json
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use json_module, only: json_file
   implicit none
   private

   public :: trc_basis_from_json

   real(dp), parameter :: PI = 3.14159265358979323846_dp

   type :: shell_t
      integer :: l = 0, np = 0
      real(dp), allocatable :: e(:), c(:)
   end type shell_t

   type :: element_t
      integer :: z = 0, nshell = 0
      type(shell_t), allocatable :: sh(:)
   end type element_t

contains

   !
   ! `at_z` atomic numbers, `at_r` (3, natm) in Bohr; the basis comes back
   ! built (host side; the caller moves it to the device).
   !
   subroutine trc_basis_from_json(path, natm, at_z, at_r, b, error)
      character(len=*), intent(in) :: path
      integer, intent(in) :: natm
      integer, intent(in) :: at_z(natm)
      real(dp), intent(in) :: at_r(3, natm)
      type(trc_basis_t), intent(out) :: b
      type(error_t), intent(inout) :: error

      type(json_file) :: json
      type(element_t), allocatable :: el(:)
      integer :: nel, ia, ie, is, nsh, maxnp, k
      integer, allocatable :: sh_l(:), sh_np(:), which(:)
      real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :), z(:)
      logical :: ok

      call json%initialize()
      call json%load(filename=path)
      if (json%failed()) then
         call error%set(ERROR_VALIDATION, "trc_basis_json: cannot read "//trim(path))
         call json%destroy()
         return
      end if

      ! One read per distinct element.
      allocate (el(natm), which(natm))
      nel = 0
      do ia = 1, natm
         ok = .false.
         do ie = 1, nel
            if (el(ie)%z == at_z(ia)) then
               which(ia) = ie; ok = .true.; exit
            end if
         end do
         if (ok) cycle
         nel = nel + 1
         which(ia) = nel
         call read_element(json, path, at_z(ia), el(nel), error)
         if (error%has_error()) then
            call json%destroy()
            return
         end if
      end do
      call json%destroy()

      nsh = 0; maxnp = 1
      do ia = 1, natm
         nsh = nsh + el(which(ia))%nshell
         do is = 1, el(which(ia))%nshell
            maxnp = max(maxnp, el(which(ia))%sh(is)%np)
         end do
      end do
      allocate (sh_l(nsh), sh_np(nsh), sh_e(maxnp, nsh), sh_c(maxnp, nsh), sh_r(3, nsh), z(natm))
      sh_e = 0.0_dp; sh_c = 0.0_dp
      k = 0
      do ia = 1, natm
         ie = which(ia)
         do is = 1, el(ie)%nshell
            k = k + 1
            sh_l(k) = el(ie)%sh(is)%l
            sh_np(k) = el(ie)%sh(is)%np
            sh_e(1:sh_np(k), k) = el(ie)%sh(is)%e
            sh_c(1:sh_np(k), k) = el(ie)%sh(is)%c
            sh_r(:, k) = at_r(:, ia)
         end do
      end do
      z = real(at_z, dp)
      call b%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, natm, z, at_r, maxnp)
   end subroutine trc_basis_from_json

   subroutine read_element(json, path, z, el, error)
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      integer, intent(in) :: z
      type(element_t), intent(out) :: el
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: key, sp, text
      integer, allocatable :: lam(:)
      integer :: nsh_json, ish, ncol, icol, ip, np, k, ios, nsh
      real(dp) :: v
      logical :: found

      el%z = z
      key = "elements."//itoa(z)
      call json%info(key//".electron_shells", found=found, n_children=nsh_json)
      if (.not. found .or. nsh_json <= 0) then
         call error%set(ERROR_VALIDATION, "trc_basis_json: element "//itoa(z)//" not in "//trim(path))
         return
      end if

      ! Count shells after splitting: one per coefficient column.
      nsh = 0
      do ish = 1, nsh_json
         call json%get(key//".electron_shells("//itoa(ish)//").angular_momentum", lam, found)
         if (.not. found) cycle
         nsh = nsh + columns(json, key, ish, size(lam))
      end do
      el%nshell = nsh
      allocate (el%sh(nsh))

      k = 0
      do ish = 1, nsh_json
         call json%get(key//".electron_shells("//itoa(ish)//").angular_momentum", lam, found)
         if (.not. found) cycle
         sp = key//".electron_shells("//itoa(ish)//")"
         call json%info(sp//".exponents", found=found, n_children=np)
         if (.not. found .or. np <= 0) then
            call error%set(ERROR_VALIDATION, "trc_basis_json: shell without exponents in "//trim(path))
            return
         end if
         ncol = columns(json, key, ish, size(lam))
         do icol = 1, ncol
            k = k + 1
            ! One l per column when the file lists one per column (SP); one
            ! l for all of them in a general contraction.
            if (size(lam) == ncol) then
               el%sh(k)%l = lam(icol)
            else
               el%sh(k)%l = lam(1)
            end if
            el%sh(k)%np = np
            allocate (el%sh(k)%e(np), el%sh(k)%c(np))
            do ip = 1, np
               call json%get(sp//".exponents("//itoa(ip)//")", text, found)
               if (.not. found) then
                  call error%set(ERROR_VALIDATION, "trc_basis_json: missing exponent in "//trim(path))
                  return
               end if
               read (text, *, iostat=ios) v
               if (ios /= 0) then
                  call error%set(ERROR_VALIDATION, "trc_basis_json: malformed exponent '"//text//"'")
                  return
               end if
               el%sh(k)%e(ip) = v
               call json%get(sp//".coefficients("//itoa(icol)//")("//itoa(ip)//")", text, found)
               if (.not. found) then
                  call error%set(ERROR_VALIDATION, "trc_basis_json: coefficients do not match exponents in "// &
                                 trim(path))
                  return
               end if
               read (text, *, iostat=ios) v
               if (ios /= 0) then
                  call error%set(ERROR_VALIDATION, "trc_basis_json: malformed coefficient '"//text//"'")
                  return
               end if
               el%sh(k)%c(ip) = v
            end do
            call normalise(el%sh(k)%l, el%sh(k)%np, el%sh(k)%e, el%sh(k)%c)
         end do
      end do
   end subroutine read_element

   ! Coefficient columns of a JSON shell, or one per l when none are listed.
   integer function columns(json, key, ish, nl) result(n)
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: key
      integer, intent(in) :: ish, nl
      logical :: found
      integer :: nc
      call json%info(key//".electron_shells("//itoa(ish)//").coefficients", found=found, n_children=nc)
      if (found .and. nc > 0) then
         n = nc
      else
         n = nl
      end if
   end function columns

   !
   ! metalquicha's normalisation (mqc_basis_normalization): each primitive
   ! to unit self-overlap, then the contraction to unit self-overlap. The
   ! coefficients come back ready for trc_basis_t%build, which folds
   ! common_fac_sp itself.
   !
   pure subroutine normalise(l, np, e, c)
      integer, intent(in) :: l, np
      real(dp), intent(in) :: e(np)
      real(dp), intent(inout) :: c(np)
      real(dp) :: dfac, ovl, s
      integer :: i, j, k
      dfac = 1.0_dp
      do k = 1, l
         dfac = dfac*real(2*k - 1, dp)
      end do
      ovl = 0.0_dp
      do i = 1, np
         do j = 1, np
            ovl = ovl + (sqrt(4.0_dp*e(i)*e(j))/(e(i) + e(j)))**(real(l, dp) + 1.5_dp)*c(i)*c(j)
         end do
      end do
      s = ovl**(-0.5_dp)
      do i = 1, np
         c(i) = s*c(i)*sqrt(2.0_dp**l/(PI**1.5_dp*dfac)*(2.0_dp*e(i))**(real(l, dp) + 1.5_dp))
      end do
   end subroutine normalise

   pure function itoa(i) result(t)
      integer, intent(in) :: i
      character(len=:), allocatable :: t
      character(len=16) :: buf
      write (buf, '(i0)') i
      t = trim(buf)
   end function itoa

end module trc_basis_json
