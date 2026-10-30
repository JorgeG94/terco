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
   use trc_api, only: trc_basis_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use json_module, only: json_file
   implicit none
   private

   public :: trc_basis_from_json, trc_fitted_shell_from_json, trc_basis_overrides_t, trc_basis_json_overrides

   !
   ! The rarely-set, never-combined-with-each-other knobs of
   ! trc_basis_from_json, bundled so the call keeps one optional dummy
   ! instead of one per knob (docs/CODE_STRUCTURE.md's arg-count limit --
   ! `force_cart` stays its own dummy since nearly every caller sets it).
   ! Build one with trc_basis_json_overrides(charges=..., prune=...), never
   ! with the structure constructor: nvfortran 26.5 leaves an omitted
   ! allocatable scalar component ALLOCATED (and .false.) when another
   ! component comes from an allocatable component of a pointer target --
   ! trc_capi's `trc_basis_overrides_t(charges=cx%at_z)` silently turned
   ! pruning off -- and aborts with an internal compiler error on a
   ! constructor mixing an allocatable and a default-initialised component.
   ! Plain assignment inside the builder is correct under both compilers.
   !
   type :: trc_basis_overrides_t
      private
      !! Private so the structure constructor is unreachable from outside
      !! this module -- the only thing that actually fixes the nvfortran
      !! defect above: a caller cannot write
      !! `trc_basis_overrides_t(charges=cx%at_z)` (or any other field-by-
      !! field construction) at all now, only `trc_basis_json_overrides`,
      !! which this module's own `trc_basis_from_json`/`read_element` read
      !! back through ordinary component access (same module, so `private`
      !! does not block it).
      real(dp), allocatable :: charges(:)
      !! Nuclear charge per atom; trc_basis_from_json's own `charges`.
      logical :: prune_set = .false.
      !! Whether `prune` overrides `TRC_PRUNE_PRIMITIVES`.
      logical :: prune = .true.
      !! The override, when `prune_set`.
   end type trc_basis_overrides_t

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
   function trc_basis_json_overrides(charges, prune) result(overrides)
      !! A trc_basis_overrides_t with just the knobs given (see the type's
      !! comment for why this, not the structure constructor).
      real(dp), intent(in), optional :: charges(:)
      logical, intent(in), optional :: prune
      type(trc_basis_overrides_t) :: overrides
      if (present(charges)) overrides%charges = charges
      if (present(prune)) then
         overrides%prune_set = .true.
         overrides%prune = prune
      end if
   end function trc_basis_json_overrides

   subroutine trc_basis_from_json(path, natm, elements, at_r, b, error, force_cart, overrides)
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
      type(trc_basis_overrides_t), intent(in), optional :: overrides
      !! The two knobs no caller sets together (trc_basis_overrides_t's own
      !! header): `charges`, nuclear charge per atom, defaulting to
      !! real(elements, dp) -- a ghost (trc_set_elements) passes 0 here for
      !! an element >= 1; and `prune`, overriding `TRC_PRUNE_PRIMITIVES`
      !! (prune_enabled()) for this one build. A caller that compares
      !! AO-level quantities (not just an energy or a gradient) against a
      !! pyscf reference built straight from this file's own literal
      !! columns -- bse_json_basis in test/pyscf_refs/generate.py, which
      !! does not remove the redundancy prune_block does -- sets
      !! `overrides%prune = .false.` to keep matching it; see that
      !! component's own call sites for which ones.

      type(json_file) :: json
      type(element_t), allocatable :: el(:)
      integer :: nel, ia, ie, is, nsh, maxnp, k
      integer, allocatable :: sh_l(:), sh_np(:), which(:)
      real(dp), allocatable :: sh_e(:, :), sh_c(:, :), sh_r(:, :), z(:)
      logical :: ok, any_sph, force, do_prune

      force = .false.
      if (present(force_cart)) force = force_cart
      do_prune = prune_enabled()
      if (present(overrides)) then
         if (overrides%prune_set) do_prune = overrides%prune
      end if

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
         call read_element(json, path, elements(ia), el(nel), error, prune=do_prune)
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
      z = real(elements, dp)
      if (present(overrides)) then
         if (allocated(overrides%charges)) z = overrides%charges
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

   subroutine read_element(json, path, z, el, error, raw, prune)
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
      logical, intent(in), optional :: prune
      !! Eliminate primitives a general contraction's column shares with a
      !! single-primitive column of the same block and angular momentum,
      !! then drop every coefficient that comes out exactly zero. Default
      !! .false. (the caller decides; `trc_basis_from_json` passes
      !! `prune_enabled()`). Never applied with `raw` -- see PRUNING, below.

      character(len=:), allocatable :: key, sp, text, ftype
      integer, allocatable :: lam(:)
      integer :: nsh_json, ish, ncol, icol, ip, np, k, ios, nsh, k0, ks
      real(dp) :: v
      logical :: found, skip_norm, found_ft, sph_ish, do_prune

      skip_norm = .false.
      if (present(raw)) skip_norm = raw
      do_prune = .false.
      if (present(prune)) do_prune = prune .and. .not. skip_norm
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
         k0 = k + 1
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
         end do
         ! PRUNING, within this one JSON block (its columns all share the
         ! same raw exponent list, which is what makes the elimination exact
         ! rather than an approximation -- see prune_block's own header).
         ! Done before `normalise` so the normalisation that follows sees
         ! each column's REAL primitive count, not the padded one.
         if (do_prune) call prune_block(el%sh(k0:k))
         do ks = k0, k
            if (.not. skip_norm) call normalise(el%sh(ks)%l, el%sh(ks)%np, el%sh(ks)%e, el%sh(ks)%c)
         end do
      end do
   end subroutine read_element

   !
   ! PRUNING a general contraction's redundant primitives.
   !
   ! cc-pVXZ (and other correlation-consistent sets) share one exponent list
   ! across a block's columns, and some of those columns are a single
   ! primitive (coefficient 1 at one exponent, 0 at the rest) that is ALSO
   ! used, with a small nonzero coefficient, inside the block's other
   ! (contracted) columns. That primitive is not independent information:
   ! the function it is part of is already in the basis set as the
   ! single-primitive column's own AO, so the contracted column's component
   ! along it can be transformed away without touching the SPAN of the
   ! block's functions -- replace contracted column A by
   ! A' = A - c_A(p)*B, B the single-primitive column with a 1 at p. B is
   ! zero everywhere else, so A' is A UNCHANGED at every other primitive and
   ! EXACTLY zero at p; no arithmetic is needed, p is just set to zero.
   !
   ! This is pyscf 2.14's own "optimized general contraction" cc-pVXZ data,
   ! reverse-engineered against the BSE files metalquicha ships: per-element
   ! per-l primitive counts after this routine match pyscf's
   ! `gto.basis.load` exactly (O/cc-pVQZ's 12-primitive s columns go to 9,
   ! matching pyscf's shipped 9; the module header this feeds,
   ! trc_decontract.F90, inherits the same win on the s/p sets of every
   ! first-row atom). w20/cc-pVQZ: 1440 primitives at <=12/shell before this
   ! routine, 1140 at <=9/shell after -- the same counts pyscf's bundled set
   ! gives for the same molecule, same energy to round-off.
   !
   ! Scoped to one JSON block (one shared exponent list) and one angular
   ! momentum within it, so an SP block's S and P columns -- which share
   ! exponents but are different operators -- never prune each other.
   !
   subroutine prune_block(sh)
      type(shell_t), intent(inout) :: sh(:)
      integer :: i, j, piv, nnz, keep
      real(dp), allocatable :: c0(:, :)
      ! A column whose nonzero primitives are EXACTLY a subset of the
      ! block's pivots (no primitive outside that set) comes out all zero
      ! below: A = sum_piv c_A(piv)*B_piv with nothing left over, so it was
      ! already linearly dependent on those pivots before this routine
      ! touched it -- the elimination is still exact. But the caller's
      ! next step, `normalise`, divides by sqrt(overlap), and an all-zero
      ! column has zero overlap, so it would come back with NaN
      ! coefficients rather than the harmless literal zero column
      ! `trc_api`'s `basis_build` already knows to leave inert. Snapshot
      ! every column here so a column that goes all zero can be restored
      ! to what BSE actually wrote instead of reaching `normalise` zeroed.
      allocate (c0(sh(1)%np, size(sh)))
      do i = 1, size(sh)
         c0(:, i) = sh(i)%c
      end do
      do i = 1, size(sh)
         ! A pivot: exactly one nonzero coefficient over the block's shared
         ! exponent list.
         nnz = count(sh(i)%c /= 0.0_dp)
         if (nnz /= 1) cycle
         piv = findloc(sh(i)%c /= 0.0_dp, .true., dim=1)
         do j = 1, size(sh)
            if (j == i) cycle
            if (sh(j)%l /= sh(i)%l) cycle
            ! Exact zero, not a subtraction: B's only nonzero is at piv with
            ! coefficient 1 there (BSE's own convention for these columns),
            ! so A - c_A(piv)*B is A with position piv set to 0 and every
            ! other position untouched by construction.
            sh(j)%c(piv) = 0.0_dp
         end do
      end do
      ! Restore any column the loop above reduced to all zeros: it is
      ! linearly dependent on the pivots it duplicated (same as it always
      ! was), not a new, empty basis function, so the conservative choice
      ! is to leave it exactly as read rather than let `normalise` turn it
      ! into NaN or silently drop a function from the span. A lindep
      ! screen downstream is what already exists to handle the resulting
      ! redundancy.
      do i = 1, size(sh)
         if (all(sh(i)%c == 0.0_dp)) sh(i)%c = c0(:, i)
      end do
      ! Compact every column of the block to its surviving primitives. A
      ! pivot column itself shrinks from the block's full exponent count
      ! down to the 1 it always was; a pruned contracted column shrinks by
      ! however many pivots it shared.
      do i = 1, size(sh)
         keep = count(sh(i)%c /= 0.0_dp)
         if (keep == sh(i)%np) cycle   ! nothing padded away
         if (keep == 0) cycle   ! defensive only: BSE itself wrote an
         ! all-zero column (the restore above only fires on one THIS
         ! routine produced) -- leave it exactly as read; `basis_build`
         ! already treats a literal all-zero shell as inert.
         block
            real(dp), allocatable :: e2(:), c2(:)
            integer :: ip, kk
            allocate (e2(keep), c2(keep))
            kk = 0
            do ip = 1, sh(i)%np
               if (sh(i)%c(ip) == 0.0_dp) cycle
               kk = kk + 1
               e2(kk) = sh(i)%e(ip)
               c2(kk) = sh(i)%c(ip)
            end do
            call move_alloc(e2, sh(i)%e)
            call move_alloc(c2, sh(i)%c)
            sh(i)%np = keep
         end block
      end do
   end subroutine prune_block

   !
   ! TRC_PRUNE_PRIMITIVES, read once per basis build. On by default; "0"
   ! turns the elimination above off, for an A/B comparison against the
   ! unpruned (BSE-literal) basis without rebuilding.
   !
   logical function prune_enabled() result(res)
      character(len=8) :: value
      value = " "
      call get_environment_variable("TRC_PRUNE_PRIMITIVES", value)
      res = (trim(value) /= "0")
   end function prune_enabled

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
      character(len=16) :: buf
      write (buf, "(i0)") i
      t = trim(buf)
   end function itoa

end module trc_basis_json
