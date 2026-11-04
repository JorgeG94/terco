!
! An effective core potential from a MolSSI Basis Set Exchange JSON file --
! the same file trc_basis_json.F90 reads `electron_shells` from, its
! `ecp_potentials`/`ecp_electrons` blocks instead.
!
! Ported from metalquicha's mqc_json_ecp_reader (mqc_ecp's shape, not its
! row layout: this module packs libcint/PySCF-style flat rows -- see
! trc_ecp_types.F90 -- so the integral driver's grouping ports unchanged).
!
! THE LOCAL CHANNEL IS THE ONE WITH THE HIGHEST l, not the one listed first
! -- BSE lists Iodine's def2-SVP potential with its l=3 (local) entry FIRST
! and the projected l=0,1,2 entries after, but nothing says every element's
! file is written that way, so this always scans and picks the maximum.
!
! A CHANNEL WHOSE PRIMITIVES MIX r-POWERS IS SPLIT INTO ONE ROW PER POWER.
! `ecp_rad_part` (trc_ecp_radial.F90) applies a row's `ecp_rpow` to every one
! of the row's primitives; a channel built from primitives of r^0 and r^2
! packed into a single row would apply r^0 (or whichever power the first
! primitive happened to carry) to all of them. def2's own ECPs here are
! uniform per channel (checked: every primitive of every def2-svp.json
! channel carries r_exponents = 2), so this split is dead code on today's
! test data and exists for tabulations that are not. An element with no
! `ecp_potentials` is not an error: it comes back with no rows for that
! atom, same as metalquicha's `has_ecp = .false.`.
!
module trc_ecp_json
   use pic_types, only: dp
   use pic_strings, only: to_string
   use trc_ecp_types, only: trc_ecp_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use json_module, only: json_file
   implicit none
   private
   public :: trc_ecp_from_json

   type :: channel_t
      integer :: l = 0
      integer :: nprim = 0
      integer, allocatable :: rpow(:)
      real(dp), allocatable :: e(:), c(:)
   end type channel_t

   type :: element_ecp_t
      integer :: z = 0
      logical :: has_ecp = .false.
      integer :: ncore = 0
      type(channel_t), allocatable :: ch(:)   !! local first, then projected
      integer :: nchannel = 0
   end type element_ecp_t

contains

   !
   ! `elements` (1..118) picks the ECP for each atom, same convention as
   ! `trc_basis_from_json`'s own `elements` argument. `ecp` comes back built
   ! (host side); an element the file has no `ecp_potentials` for leaves
   ! that atom with `at_ncore = 0` and no rows, not an error.
   !
   subroutine trc_ecp_from_json(path, natm, elements, ecp, error)
      character(len=*), intent(in) :: path
      integer, intent(in) :: natm
      integer, intent(in) :: elements(natm)
      type(trc_ecp_t), intent(out) :: ecp
      type(error_t), intent(inout) :: error

      type(json_file) :: json
      type(element_ecp_t), allocatable :: el(:)
      integer, allocatable :: which(:)
      integer :: nel, ia, ie, nrow, maxnp, k, ic

      call json%initialize()
      call json%load(filename=path)
      if (json%failed()) then
         call error%set(ERROR_VALIDATION, "trc_ecp_json: cannot read "//trim(path))
         call json%destroy()
         return
      end if

      allocate (el(natm), which(natm))
      nel = 0
      do ia = 1, natm
         ie = find_element(el, nel, elements(ia))
         if (ie == 0) then
            nel = nel + 1
            ie = nel
            call read_element_ecp(json, path, elements(ia), el(ie), error)
            if (error%has_error()) then
               call json%destroy()
               return
            end if
         end if
         which(ia) = ie
      end do
      call json%destroy()

      nrow = 0
      maxnp = 1
      do ia = 1, natm
         do ic = 1, el(which(ia))%nchannel
            nrow = nrow + 1
            maxnp = max(maxnp, el(which(ia))%ch(ic)%nprim)
         end do
      end do

      ecp%nrow = nrow
      ecp%natm = natm
      ecp%maxnp = maxnp
      allocate (ecp%ecp_at(nrow), ecp%ecp_l(nrow), ecp%ecp_rpow(nrow), ecp%ecp_np(nrow))
      allocate (ecp%ecp_e(maxnp, nrow), ecp%ecp_c(maxnp, nrow))
      allocate (ecp%at_ncore(natm))
      ecp%ecp_e = 0.0_dp
      ecp%ecp_c = 0.0_dp
      ecp%at_ncore = 0
      ecp%has_ecp = .false.

      k = 0
      do ia = 1, natm
         ie = which(ia)
         if (.not. el(ie)%has_ecp) cycle
         ecp%at_ncore(ia) = el(ie)%ncore
         ecp%has_ecp = .true.
         do ic = 1, el(ie)%nchannel
            k = k + 1
            ecp%ecp_at(k) = ia
            ecp%ecp_l(k) = el(ie)%ch(ic)%l
            ! Every primitive of a row shares one r-power by construction
            ! (split_by_rpow built it that way); rpow(1) speaks for the row.
            ecp%ecp_rpow(k) = el(ie)%ch(ic)%rpow(1)
            ecp%ecp_np(k) = el(ie)%ch(ic)%nprim
            ecp%ecp_e(1:el(ie)%ch(ic)%nprim, k) = el(ie)%ch(ic)%e
            ecp%ecp_c(1:el(ie)%ch(ic)%nprim, k) = el(ie)%ch(ic)%c
         end do
      end do
   end subroutine trc_ecp_from_json

   integer function find_element(el, nel, z) result(idx)
      type(element_ecp_t), intent(in) :: el(:)
      integer, intent(in) :: nel, z
      integer :: i
      idx = 0
      do i = 1, nel
         if (el(i)%z == z) then
            idx = i
            return
         end if
      end do
   end function find_element

   subroutine read_element_ecp(json, path, z, el, error)
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: path
      integer, intent(in) :: z
      type(element_ecp_t), intent(out) :: el
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: key, ptype
      integer, allocatable :: chan_l(:)
      type(channel_t), allocatable :: raw(:), split(:)
      integer :: npot, ipot, local_index, nchan, ic, nsplit, ntot, off
      logical :: found

      el%z = z
      key = "elements."//to_string(z)
      call json%info(key//".ecp_potentials", found=found, n_children=npot)
      if (.not. found .or. npot <= 0) return   ! no ECP for this element

      call json%get(key//".ecp_electrons", el%ncore, found)
      if (.not. found) then
         call error%set(ERROR_VALIDATION, "trc_ecp_json: element "//to_string(z)// &
                        " has ecp_potentials but no ecp_electrons in "//trim(path))
         return
      end if

      allocate (chan_l(npot), raw(npot))
      do ipot = 1, npot
         call json%get(key//".ecp_potentials("//to_string(ipot)//").ecp_type", ptype, found)
         if (found) then
            if (trim(ptype) /= "scalar_ecp") then
               call error%set(ERROR_VALIDATION, "trc_ecp_json: element "//to_string(z)// &
                              " carries a '"//trim(ptype)//"' ECP channel in "//trim(path)// &
                              "; only scalar_ecp is supported (no spin-orbit)")
               return
            end if
         end if
         call read_raw_channel(json, key, ipot, raw(ipot), path, error)
         if (error%has_error()) return
         chan_l(ipot) = raw(ipot)%l
      end do
      local_index = maxloc(chan_l, dim=1)

      ! Local channel first, then every projected one, each split into one
      ! row per distinct r-power its primitives carry.
      call split_by_rpow(raw(local_index), split, nsplit)
      ntot = nsplit
      do ipot = 1, npot
         if (ipot == local_index) cycle
         call split_by_rpow(raw(ipot), split, nsplit)
         ntot = ntot + nsplit
      end do

      allocate (el%ch(ntot))
      off = 0
      call split_by_rpow(raw(local_index), split, nsplit)
      do ic = 1, nsplit
         el%ch(off + ic) = split(ic)
         el%ch(off + ic)%l = -1   ! the local-channel sentinel, pyscf's own
      end do
      off = off + nsplit
      do ipot = 1, npot
         if (ipot == local_index) cycle
         call split_by_rpow(raw(ipot), split, nsplit)
         el%ch(off + 1:off + nsplit) = split(1:nsplit)
         off = off + nsplit
      end do
      el%nchannel = ntot
      el%has_ecp = .true.
   end subroutine read_element_ecp

   !
   ! One ecp_potentials() entry, exactly as the file lists it: one l, and
   ! parallel r_exponents/gaussian_exponents/coefficients(1) lists.
   !
   subroutine read_raw_channel(json, key, ipot, ch, path, error)
      type(json_file), intent(inout) :: json
      character(len=*), intent(in) :: key
      integer, intent(in) :: ipot
      type(channel_t), intent(out) :: ch
      character(len=*), intent(in) :: path
      type(error_t), intent(inout) :: error

      character(len=:), allocatable :: cpath, text
      integer, allocatable :: momenta(:), powers(:)
      integer :: nprim, iprim, ios, nc
      real(dp) :: v
      logical :: found

      cpath = key//".ecp_potentials("//to_string(ipot)//")"
      call json%get(cpath//".angular_momentum", momenta, found)
      if (found .and. allocated(momenta)) then
         if (size(momenta) < 1) found = .false.
      else
         found = .false.
      end if
      if (.not. found) then
         call error%set(ERROR_VALIDATION, "trc_ecp_json: ECP channel without "// &
                        "angular_momentum in "//trim(path))
         return
      end if
      ch%l = momenta(1)

      call json%get(cpath//".r_exponents", powers, found)
      if (.not. found .or. .not. allocated(powers)) then
         call error%set(ERROR_VALIDATION, "trc_ecp_json: ECP channel without "// &
                        "r_exponents in "//trim(path))
         return
      end if
      nprim = size(powers)

      call json%info(cpath//".gaussian_exponents", found=found, n_children=nc)
      if (.not. found .or. nc /= nprim) then
         call error%set(ERROR_VALIDATION, "trc_ecp_json: ECP channel with mismatched "// &
                        "exponent count in "//trim(path))
         return
      end if

      ch%nprim = nprim
      allocate (ch%rpow(nprim), ch%e(nprim), ch%c(nprim))
      ch%rpow = powers
      do iprim = 1, nprim
         call json%get(cpath//".gaussian_exponents("//to_string(iprim)//")", text, found)
         if (.not. found) then
            call error%set(ERROR_VALIDATION, "trc_ecp_json: missing ECP exponent in "//trim(path))
            return
         end if
         read (text, *, iostat=ios) v
         if (ios /= 0) then
            call error%set(ERROR_VALIDATION, "trc_ecp_json: malformed ECP exponent '"// &
                           text//"' in "//trim(path))
            return
         end if
         ch%e(iprim) = v

         call json%get(cpath//".coefficients(1)("//to_string(iprim)//")", text, found)
         if (.not. found) then
            call error%set(ERROR_VALIDATION, "trc_ecp_json: missing ECP coefficient in "//trim(path))
            return
         end if
         read (text, *, iostat=ios) v
         if (ios /= 0) then
            call error%set(ERROR_VALIDATION, "trc_ecp_json: malformed ECP coefficient '"// &
                           text//"' in "//trim(path))
            return
         end if
         ch%c(iprim) = v
      end do
   end subroutine read_raw_channel

   !
   ! Split one channel's primitives into one row per DISTINCT r-power,
   ! preserving the order those powers first appear in. For every channel
   ! in today's def2-svp.json this produces exactly one row (see this
   ! module's header); it is still applied unconditionally so a tabulation
   ! that does mix powers is handled rather than silently wrong.
   !
   subroutine split_by_rpow(ch, out, nout)
      type(channel_t), intent(in) :: ch
      type(channel_t), allocatable, intent(out) :: out(:)
      integer, intent(out) :: nout

      integer, allocatable :: distinct(:), group(:)
      integer :: i, j, g, n

      allocate (distinct(ch%nprim), group(ch%nprim))
      nout = 0
      do i = 1, ch%nprim
         g = 0
         do j = 1, nout
            if (distinct(j) == ch%rpow(i)) then
               g = j
               exit
            end if
         end do
         if (g == 0) then
            nout = nout + 1
            distinct(nout) = ch%rpow(i)
            g = nout
         end if
         group(i) = g
      end do

      allocate (out(nout))
      do g = 1, nout
         n = count(group == g)
         out(g)%l = ch%l
         out(g)%nprim = n
         allocate (out(g)%rpow(n), out(g)%e(n), out(g)%c(n))
         out(g)%rpow = distinct(g)
         n = 0
         do i = 1, ch%nprim
            if (group(i) /= g) cycle
            n = n + 1
            out(g)%e(n) = ch%e(i)
            out(g)%c(n) = ch%c(i)
         end do
      end do
      deallocate (distinct, group)
   end subroutine split_by_rpow

end module trc_ecp_json
