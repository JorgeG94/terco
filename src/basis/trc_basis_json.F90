!
! A basis from a MolSSI Basis Set Exchange JSON file, the files metalquicha
! ships under basis_sets/, so a terco run can take the orbital and auxiliary
! sets mqc would hand it without mqc in the process.
!
! Ported from metalquicha's mqc_json_basis_reader
!
module trc_basis_json
   use pic_types, only: dp
   use trc_constants, only: PI
   use trc_basis_types, only: trc_basis_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use json_module, only: json_file
   implicit none
   private

   public :: trc_basis_from_json, trc_fitted_shell_from_json

   type :: shell_t
      integer :: l = 0, np = 0
      real(dp), allocatable :: e(:), c(:)
      logical :: sph = .false.   !! this shell's own function_type said "gto_spherical"
   end type shell_t

   type :: element_t
      integer :: z = 0, nshell = 0
      type(shell_t), allocatable :: sh(:)
   end type element_t

contains

   !
   ! `elements` (1..118) picks the basis for each atom; `at_r` (3, natm) is
   ! in Bohr. `charges` is the nuclear charge the built basis carries as
   ! `at_z`, defaulting to `elements` itself (every caller before ghost
   ! atoms existed): a ghost passes an element with charge 0, so its shells
   ! are looked up and placed exactly as a real atom of that element's, but
   ! it contributes no nuclear charge and (through `at_z`) no electron.
   ! The basis comes back built (host side; the caller moves it to the
   ! device).
   !
   subroutine trc_basis_from_json(path, natm, elements, at_r, b, error, force_cart, charges)
      character(len=*), intent(in) :: path
      integer, intent(in) :: natm
      integer, intent(in) :: elements(natm)
      real(dp), intent(in) :: at_r(3, natm)
      type(trc_basis_t), intent(out) :: b
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: force_cart
      !! Force the Cartesian AO space regardless of what the file says --
      !! trc_scf_options_t%force_cartesian's route in. Default .false.:
      !! follow the basis, which is what every existing caller gets.
      real(dp), intent(in), optional :: charges(natm)
      !! Nuclear charge per atom; defaults to real(elements, dp). A ghost
      !! (trc_set_elements) passes 0 here for an element >= 1.

      type(json_file) :: json
      type(element_t), allocatable :: el(:)
      integer :: nel, ia, ie, is, nsh, maxnp, k
      integer, allocatable :: sh_l(:), sh_np(:), which(:)
      real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :), z(:)
      logical :: ok, any_sph, force

      force = .false.
      if (present(force_cart)) force = force_cart

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
            if (el(ie)%z == elements(ia)) then
               which(ia) = ie
               ok = .true.
               exit
            end if
         end do
         if (ok) cycle
         nel = nel + 1
         which(ia) = nel
         call read_element(json, path, elements(ia), el(nel), error)
         if (error%has_error()) then
            call json%destroy()
            return
         end if
      end do
      call json%destroy()

      nsh = 0
      maxnp = 1
      do ia = 1, natm
         nsh = nsh + el(which(ia))%nshell
         do is = 1, el(which(ia))%nshell
            maxnp = max(maxnp, el(which(ia))%sh(is)%np)
         end do
      end do
      allocate (sh_l(nsh), sh_np(nsh), sh_e(maxnp, nsh), sh_c(maxnp, nsh), sh_r(3, nsh), z(natm))
      sh_e = 0.0_dp
      sh_c = 0.0_dp
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
      if (present(charges)) then
         z = charges
      else
         z = real(elements, dp)
      end if
      call b%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, natm, z, at_r, maxnp, at_elem=elements)

      ! Any shell of any element tagged "gto_spherical" is enough: a real
      ! basis set is one family or the other throughout its l>=2 shells
      ! (see this module's own header), so one true tag settles the whole
      ! basis rather than needing every shell to agree.
      any_sph = .false.
      do ie = 1, nel
         do is = 1, el(ie)%nshell
            if (el(ie)%sh(is)%sph) any_sph = .true.
         end do
      end do
      b%spherical = any_sph .and. .not. force
   end subroutine trc_basis_from_json

   !
   ! One element's fitted-potential shell, exactly as the file lists it.
   !
   ! For the SAP sets on the BSE (sap_helfem_large and friends, role
   ! "guess"). Their coefficients are the CHARGES of a fitted distribution,
   ! not wavefunction coefficients: they sum to -Z for every element, which
   ! is the check that says so and the one this routine makes. Normalising
   ! them to unit overlap, which is what every other basis in this file
   ! wants, would be wrong in a way that leaves the numbers looking sane.
   !
   subroutine trc_fitted_shell_from_json(path, z, e, c, error)
      character(len=*), intent(in) :: path
      integer, intent(in) :: z
      real(dp), allocatable, intent(out) :: e(:), c(:)
      type(error_t), intent(inout) :: error

      type(json_file) :: json
      type(element_t) :: el
      real(dp) :: qsum

      call json%initialize()
      call json%load_file(filename=trim(path))
      if (json%failed()) then
         call error%set(ERROR_VALIDATION, "trc_basis_json: cannot read "//trim(path))
         call json%destroy()
         return
      end if
      call read_element(json, path, z, el, error, raw=.true.)
      call json%destroy()
      if (error%has_error()) return
      if (el%nshell /= 1 .or. el%sh(1)%l /= 0) then
         call error%set(ERROR_VALIDATION, "trc_basis_json: "//trim(path)//" gives element "// &
                        itoa(z)//" more than one shell, or a shell that is not s; a fitted "// &
                        "potential is one spherical charge distribution per atom")
         return
      end if
      e = el%sh(1)%e
      c = el%sh(1)%c
      ! The coefficients are charges and must account for every electron. A
      ! set that fails this is not a fitted potential, whatever it is called,
      ! and going on would give a silently weak or wrong-signed guess.
      qsum = sum(c)
      if (abs(qsum + real(z, dp)) > 1.0e-6_dp*max(real(z, dp), 1.0_dp)) then
         call error%set(ERROR_VALIDATION, "trc_basis_json: the fitted charges for element "// &
                        itoa(z)//" in "//trim(path)//" sum to "//rtoa(qsum)//", not -"//itoa(z))
         return
      end if
   end subroutine trc_fitted_shell_from_json

   subroutine read_element(json, path, z, el, error, raw)
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      integer, intent(in) :: z
      type(element_t), intent(out) :: el
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: raw
      !! Skip the wavefunction normalisation. For a set whose coefficients
      !! are not wavefunction coefficients -- a fitted potential lists
      !! CHARGES, and normalising those to unit overlap would be a category
      !! error that still produces plausible numbers.

      character(len=:), allocatable :: key, sp, text, ftype
      integer, allocatable :: lam(:)
      integer :: nsh_json, ish, ncol, icol, ip, np, k, ios, nsh
      real(dp) :: v
      logical :: found, skip_norm, found_ft, sph_ish

      skip_norm = .false.
      if (present(raw)) skip_norm = raw
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
         ! BSE's own per-shell tag ("gto", "gto_spherical" or
         ! "gto_cartesian"): read once per JSON shell and carried onto
         ! every column it splits into. Missing entirely (an older or
         ! hand-written file) reads as Cartesian, same as "gto" or
         ! "gto_cartesian" -- today's behaviour for a file that never had
         ! an opinion.
         call json%get(sp//".function_type", ftype, found_ft)
         sph_ish = found_ft
         if (sph_ish) sph_ish = (trim(ftype) == "gto_spherical")
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
            el%sh(k)%sph = sph_ish
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
            if (.not. skip_norm) call normalise(el%sh(k)%l, el%sh(k)%np, el%sh(k)%e, el%sh(k)%c)
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
   ! to unit self-overlap, then the contraction to unit self-overlap. That
   ! is BSE's own convention (also gmshpc's and pyscf's): it makes <s|s> = 1
   ! for a fully-built shell.
   !
   ! terco's convention is not that one. `trc_basis_t%build` unconditionally
   ! folds libcint's `common_fac_sp` into every s and p shell's coefficients
   ! (trc_api.F90, duplicated in trc_ao_norm.F90 and trc_batch.F90 for the
   ! same reason it is duplicated here rather than exported: it is three
   ! lines, and the alternative widens trc_api's public surface for one
   ! boundary conversion), which makes terco's OWN raw <s|s> come out
   ! 1/(4*pi), not 1. Everything else that ever builds a shell -- the
   ! hardcoded 6-31G in test/trc_test_basis.F90, trc_decontract.F90's split
   ! shells -- hands `build` a coefficient with NO common_fac_sp in it, by
   ! construction, so the single fold inside `build` lands on the right
   ! number.
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
         c(i) = s*c(i)*sqrt(2.0_dp**l/(PI**1.5_dp*dfac)*(2.0_dp*e(i))**(real(l, dp) + 1.5_dp)) &
                /common_fac_sp(l)
      end do
   end subroutine normalise

   pure real(dp) function common_fac_sp(l) result(res)
      !! Duplicated from trc_api.F90 (and trc_ao_norm.F90, and trc_batch.F90)
      !! rather than imported -- see those modules' own comments on the same
      !! duplication. `build` folds this in once; this routine has to undo the
      !! complete BSE/pyscf normalisation's OWN copy of it first, so that the
      !! two folds net to exactly one.
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

   pure function rtoa(x) result(t)
      real(dp), intent(in) :: x
      character(len=24) :: t
      write (t, "(f0.6)") x
   end function rtoa

   pure function itoa(i) result(t)
      integer, intent(in) :: i
      character(len=:), allocatable :: t
      character(len=16) :: buffer
      write (buffer, "(i0)") i
      t = trim(buffer)
   end function itoa

end module trc_basis_json
