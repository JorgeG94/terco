!
! The effective-core-potential data model: one flat set of ECP "rows",
! structure of arrays, kept entirely separate from `trc_basis_t`'s orbital
! shells.
!
! Each row is one (atom, l, r-power) GROUP of primitives, the same packing
! PySCF/libcint use internally (`ecpbas`): `ecp_l == -1` marks the LOCAL
! channel (the one with the highest angular momentum in the potential, not
! necessarily the first one listed); every other row is a PROJECTED channel
! acting through its own l. A channel whose primitives carry more than one
! r-power is split into one row per power -- see trc_ecp_json.F90's reader
! for why that split is not optional: a single row's radial routine
! (trc_ecp_radial.F90's `ecp_rad_part`) applies ONE power to every primitive
! in the row.
!
! Rows belonging to the same (atom, l) are written consecutively by the
! reader (one channel's rows, in sequence, before the next channel's), which
! is what lets the integral driver (src/eri/trc_ecp.F90) find a channel's
! whole primitive set as one contiguous row range.
!
module trc_ecp_types
   use pic_types, only: dp
   implicit none
   private
   public :: trc_ecp_t, trc_ecp_subset_atom

   type :: trc_ecp_t
      integer :: nrow = 0, natm = 0, maxnp = 0
      logical :: has_ecp = .false.
      integer, allocatable :: ecp_at(:)     !! (nrow) atom index, 1-based
      integer, allocatable :: ecp_l(:)      !! (nrow) channel l, or -1 for local
      integer, allocatable :: ecp_rpow(:)   !! (nrow) r exponent, ONE per row
      integer, allocatable :: ecp_np(:)     !! (nrow) primitives in this row
      real(dp), allocatable :: ecp_e(:, :)  !! (maxnp, nrow) exponents
      real(dp), allocatable :: ecp_c(:, :)  !! (maxnp, nrow) coefficients, raw
      !! -- an ECP channel is not an orbital shell and carries no
      !! common_fac_sp; the orbital shell's own factor (already folded into
      !! trc_basis_t%sh_c once, in basis construction) is what the driver
      !! multiplies in, so this is read verbatim from the BSE JSON.
      integer, allocatable :: at_ncore(:)   !! (natm) core electrons replaced,
      !! 0 for an atom with no ECP (phase 3's electron count and nuclear
      !! charge contract read this; nothing in phases 1-2 does).
   contains
      procedure :: release => ecp_release
   end type trc_ecp_t

contains

   subroutine ecp_release(this)
      class(trc_ecp_t), intent(inout) :: this
      if (allocated(this%ecp_at)) deallocate (this%ecp_at)
      if (allocated(this%ecp_l)) deallocate (this%ecp_l)
      if (allocated(this%ecp_rpow)) deallocate (this%ecp_rpow)
      if (allocated(this%ecp_np)) deallocate (this%ecp_np)
      if (allocated(this%ecp_e)) deallocate (this%ecp_e)
      if (allocated(this%ecp_c)) deallocate (this%ecp_c)
      if (allocated(this%at_ncore)) deallocate (this%at_ncore)
      this%nrow = 0
      this%natm = 0
      this%maxnp = 0
      this%has_ecp = .false.
   end subroutine ecp_release

   !
   ! The rows belonging to ONE atom of `this`, renumbered onto a one-atom
   ! ECP (`out%ecp_at` all 1, `out%natm = 1`) -- the same "subset" idea
   ! `trc_basis_t%subset` already uses for SAD's free-atom sub-SCF, so that
   ! sub-SCF's own H_core can carry THIS atom's U_ECP instead of running
   ! valence-only-V_ne with no core potential at all (trc_sad.F90's own
   ! header: a free atom's electron count is already right from `b%at_z`,
   ! but with no U its orbitals come out far too diffuse). `out` comes back
   ! `has_ecp = .false.` (its default) when `ia` carries no ECP rows,
   ! including every atom of an ECP-free run.
   !
   pure subroutine trc_ecp_subset_atom(this, ia, out)
      type(trc_ecp_t), intent(in) :: this
      integer, intent(in) :: ia
      type(trc_ecp_t), intent(out) :: out
      integer :: i, k, nrow

      if (.not. this%has_ecp .or. this%nrow == 0) return
      nrow = count(this%ecp_at(1:this%nrow) == ia)
      if (nrow == 0) return

      out%nrow = nrow
      out%natm = 1
      out%maxnp = this%maxnp
      allocate (out%ecp_at(nrow), out%ecp_l(nrow), out%ecp_rpow(nrow), out%ecp_np(nrow))
      allocate (out%ecp_e(this%maxnp, nrow), out%ecp_c(this%maxnp, nrow))
      allocate (out%at_ncore(1))
      k = 0
      do i = 1, this%nrow
         if (this%ecp_at(i) /= ia) cycle
         k = k + 1
         out%ecp_at(k) = 1
         out%ecp_l(k) = this%ecp_l(i)
         out%ecp_rpow(k) = this%ecp_rpow(i)
         out%ecp_np(k) = this%ecp_np(i)
         out%ecp_e(:, k) = this%ecp_e(:, i)
         out%ecp_c(:, k) = this%ecp_c(:, i)
      end do
      out%at_ncore(1) = this%at_ncore(ia)
      out%has_ecp = .true.
   end subroutine trc_ecp_subset_atom

end module trc_ecp_types
