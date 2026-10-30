!
! The basis and the screened shell-pair list: the two containers every
! integral path reads.
!
! `trc_basis_t` is a structure of arrays built once per geometry, from flat
! per-shell arrays (`build`) or libcint's `atm`/`bas`/`env` layout
! (`from_libcint`), with libcint's common_fac_sp folded into the stored
! coefficients exactly once, here. `trc_pairlist_t` is the shared screened
! shell-pair list built from a basis and a threshold, with its primitive-pair
! data precomputed once per pair. Both move to the device with `to_device`
! and leave it with `release`.
!
! A core module, below everything that consumes the two types: the basis
! reader, the integral wrappers, the Fock builds, XC, PCM, RI-MP2 and the SCF.
! `trc_api` re-exports both, so callers of the public API see no change.
!
! `build` pushes the Boys table to the device (boys_init) as it always has.
! The Rys root tables are not pushed here, since they belong to the eri
! layer: every launcher of a kernel that reads them pushes them itself
! (trc_integrals' DF wrappers, trc_eri's build, trc_scf_grad_df).
!
module trc_basis_types
   use pic_types, only: dp
   use trc_constants, only: PI
   use trc_boys, only: boys_init
   implicit none
   private
   public :: trc_basis_t, trc_pairlist_t

   ! libcint's packed layout, so a caller can pass its arrays straight in.
   integer, parameter :: ATM_SLOTS = 6, BAS_SLOTS = 8
   integer, parameter :: CHARGE_OF = 0, PTR_COORD = 1
   integer, parameter :: ATOM_OF = 0, ANG_OF = 1, NPRIM_OF = 2
   integer, parameter :: PTR_EXP = 5, PTR_COEFF = 6

   !
   ! Structure of arrays, not array of structures: the kernels read one field
   ! across many shells at once, so the fields are what should be contiguous.
   !
   type :: trc_basis_t
      integer :: nshell = 0, nao = 0, natm = 0, maxnp = 0
      integer,  allocatable :: sh_l(:)      !! angular momentum
      integer,  allocatable :: sh_np(:)     !! primitives
      integer,  allocatable :: sh_ao(:)     !! first AO of the shell, 1-based
      real(dp), allocatable :: sh_e(:, :)   !! (maxnp, nshell) exponents
      real(dp), allocatable :: sh_c(:, :)   !! (maxnp, nshell) coefficients,
                                            !! with common_fac_sp folded in
      real(dp), allocatable :: sh_r(:, :)   !! (3, nshell) centres
      integer,  allocatable :: sh_at(:)     !! (nshell) atom index of the shell's
                                            !! centre, for the gradient's shell ->
                                            !! atom map; built in `build` by an
                                            !! exact match of sh_r against at_r
      real(dp), allocatable :: at_z(:)      !! (natm) charges -- 0 for a ghost
      real(dp), allocatable :: at_r(:, :)   !! (3, natm) positions
      integer,  allocatable :: at_elem(:)   !! (natm) element (1..118), for basis
                                            !! lookup, grid sizes and guess shape;
                                            !! equal to nint(at_z) unless a caller
                                            !! placed a ghost (trc_set_elements)
      !
      ! THE SPHERICAL AO SPACE. Computed unconditionally -- it is a cheap
      ! by-product of the same loop that sets `sh_ao`, and a caller that
      ! never asks for it pays nothing extra for having it sit there.
      ! `spherical` is the only field that actually changes anything: it is
      ! what trc_scf_driver reads to decide whether the SCF runs its own
      ! matrices (S, H, D, F, C, eps) at `nao` or at `nao_sph`. The
      ! integrals themselves never see this flag -- every kernel, generator
      ! and screen stays Cartesian, exactly as today, and `nao` keeps
      ! meaning "Cartesian functions" everywhere it already did.
      !
      ! `sh_ao_sph` is `sh_ao`'s own first-AO-of-the-shell offset, but into
      ! the SPHERICAL numbering (2l+1 wide per shell instead of
      ! (l+1)(l+2)/2): trc_sph.F90's block-sparse transforms are the only
      ! reader, and need both offsets side by side to copy one shell's
      ! block at a time.
      !
      logical :: spherical = .false.
      integer :: nao_sph = 0
      integer,  allocatable :: sh_ao_sph(:)
      logical :: on_device = .false.
   contains
      !! The primary constructor: flat arrays, one entry per shell. A caller
      !! that has a basis set has these; only a caller that already links
      !! libcint has `atm`/`bas`/`env`.
      procedure :: build        => basis_build
      procedure :: subset       => basis_subset
      procedure :: from_libcint => basis_from_libcint
      procedure :: to_device    => basis_to_device
      procedure :: release      => basis_release
   end type trc_basis_t

   !
   ! Screened shell-pair list.
   !
   ! The four-centre path has had one of these since early on (`pair_bins_t`,
   ! binned by (la, lb, K, size class) with a Schwarz pre-screen); the
   ! one-electron and density-fitting paths had nothing and walked the full
   ! nshell^2 square. This is the shared object: built once per geometry from
   ! a basis and a threshold, consumed by every operator.
   !
   ! It carries TWO bounds, because Schwarz is the right bound for
   ! two-electron work and the wrong one for overlap. Both are rigorous upper
   ! bounds on the contracted block, summed over primitive pairs:
   !
   !     sbound = sum_ij |c_i c_j| (pi/zeta_ij)^(3/2) exp(-xi_ij R_AB^2)
   !     vbound = sum_ij |c_i c_j| (2 pi/zeta_ij)     exp(-xi_ij R_AB^2)
   !
   ! sbound bounds |S|; vbound bounds |V| once multiplied by the total nuclear
   ! charge, since the Boys function never exceeds 1.
   !
   type :: trc_pairlist_t
      integer :: npair = 0, nshell = 0
      real(dp) :: thresh = 0.0_dp
      integer,  allocatable :: p_i(:), p_j(:)   !! shell indices, 1-based
      real(dp), allocatable :: p_s(:), p_v(:)   !! the two bounds
      integer,  allocatable :: p_n(:)           !! primitive pairs in each
      !
      ! Primitive-pair data, computed ONCE per shell pair.
      !
      ! This is what makes the three-centre path affordable. Its kernel is
      ! called once per (shell pair, auxiliary shell), and it was rebuilding
      ! zeta, P and K_ab from the exponents on every one of those calls --
      ! so the bra product was recomputed once per auxiliary SHELL, several
      ! hundred times over for a real fitting basis. Precomputed here, it is
      ! read.
      !
      integer :: npp = 0
      integer,  allocatable :: pp_off(:)     !! (npair) first primitive pair - 1
      integer,  allocatable :: pp_n(:)       !! (npair) how many
      real(dp), allocatable :: pp_z(:)       !! (npp) zeta = alpha + beta
      real(dp), allocatable :: pp_k(:)       !! (npp) c_i c_j exp(-xi R_AB^2)
      real(dp), allocatable :: pp_p(:, :)    !! (3, npp) P
      real(dp), allocatable :: pp_a(:, :)    !! (3, npp) P - A
      logical :: on_device = .false.
   contains
      procedure :: build     => pairs_build
      procedure :: to_device => pairs_to_device
      procedure :: release   => pairs_release
   end type trc_pairlist_t

contains

   pure integer function ncart(l) result(res)
      integer, intent(in) :: l
      res = (l + 1)*(l + 2)/2
   end function ncart

   !
   ! libcint applies this per shell OUTSIDE `env`, so a container built from
   ! `env` alone is wrong by it. Folding it into sh_c here means the kernels
   ! never have to know.
   !
   pure real(dp) function common_fac_sp(l) result(res)
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

   !
   ! Build from flat per-shell arrays.
   !
   ! `normalised` says whether the caller's coefficients already carry the
   ! primitive normalisation. It is an explicit argument rather than an
   ! assumption because getting it wrong makes every integral wrong by a
   ! smooth factor that still looks plausible -- which has happened twice in
   ! this repository already.
   !
   subroutine basis_build(this, nshell, sh_l, sh_np, sh_e, sh_c, sh_r, &
                          natm, at_z, at_r, maxnp, at_elem)
      class(trc_basis_t), intent(inout) :: this
      integer,  intent(in) :: nshell, natm, maxnp
      integer,  intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in) :: sh_r(3, nshell)
      real(dp), intent(in) :: at_z(natm), at_r(3, natm)
      integer,  intent(in), optional :: at_elem(natm)
      !! Element per atom, for a ghost (at_z == 0, this build's own element
      !! still >= 1) -- absent, every caller before ghosts existed, this
      !! defaults to nint(at_z), which IS the element for a real atom.

      integer :: i, k, n, ia

      call this%release()
      call boys_init()

      this%nshell = nshell
      this%natm   = natm
      this%maxnp  = maxnp

      allocate (this%sh_l(nshell), this%sh_np(nshell), this%sh_ao(nshell))
      allocate (this%sh_ao_sph(nshell))
      allocate (this%sh_e(maxnp, nshell), this%sh_c(maxnp, nshell))
      allocate (this%sh_r(3, nshell), this%sh_at(nshell))
      allocate (this%at_z(natm), this%at_r(3, natm), this%at_elem(natm))
      this%sh_e = 0.0_dp
      this%sh_c = 0.0_dp

      this%at_z = at_z
      this%at_r = at_r
      if (present(at_elem)) then
         this%at_elem = at_elem
      else
         this%at_elem = nint(at_z)
      end if

      this%nao = 0
      this%nao_sph = 0
      this%spherical = .false.
      do i = 1, nshell
         this%sh_l(i)  = sh_l(i)
         this%sh_np(i) = sh_np(i)
         this%sh_r(:, i) = sh_r(:, i)
         ! Every caller sets sh_r(:,i) to a bit-identical copy of some
         ! at_r(:,ia) -- trc_basis_json.F90, trc_decontract.F90 (a split
         ! shell keeps its parent's centre) and basis_from_libcint all copy
         ! the coordinate, never recompute it -- so an exact match is a
         ! robust shell -> atom map, not a tolerance-guarded heuristic.
         ! sh_at stays 0 for a shell on no atom. That is legitimate here --
         ! a ghost or probe function (sapcmp's far-field check) needs no
         ! atom for the one-electron or Fock paths -- so it is
         ! eri%gradient, the one consumer of sh_at, that refuses it.
         this%sh_at(i) = 0
         do ia = 1, natm
            if (sh_r(1, i) == at_r(1, ia) .and. sh_r(2, i) == at_r(2, ia) &
                .and. sh_r(3, i) == at_r(3, ia)) then
               this%sh_at(i) = ia
               exit
            end if
         end do
         ! Primitives with a zero coefficient are dropped here, once. A
         ! general-contraction table lists every exponent of the block under
         ! every column, so the uncontracted column of a cc-pVDZ carbon s
         ! arrives with eight dead primitives out of nine; kept, they cost
         ! the quartet loops their full price for nothing -- four times the
         ! (ss|ss) work on that basis.
         n = 0
         do k = 1, sh_np(i)
            if (sh_c(k, i) == 0.0_dp) cycle
            n = n + 1
            this%sh_e(n, i) = sh_e(k, i)
            this%sh_c(n, i) = sh_c(k, i)*common_fac_sp(sh_l(i))
         end do
         if (n == 0) then
            n = sh_np(i)
            this%sh_e(1:n, i) = sh_e(1:n, i)
            this%sh_c(1:n, i) = 0.0_dp
         end if
         this%sh_np(i) = n
         this%sh_ao(i) = this%nao + 1
         this%nao = this%nao + ncart(sh_l(i))
         this%sh_ao_sph(i) = this%nao_sph + 1
         this%nao_sph = this%nao_sph + 2*sh_l(i) + 1
      end do
   end subroutine basis_build

   !
   ! Convenience adapter for callers that already hold libcint's layout.
   ! NOT the primary interface -- it unpacks into raw arrays and delegates to
   ! `build`, so there is one implementation of the import rules (notably the
   ! common_fac_sp folding) rather than two that can drift.
   !
   !
   ! Shells s0..s1 of this basis as a basis of their own, functions
   ! renumbered from 1, atoms kept. Copied field by field: the stored
   ! coefficients already carry common_fac_sp, and rebuilding from them
   ! through `build` would fold it in a second time -- which is exactly
   ! what happened to the first RI-MP2 auxiliary block, and the third time
   ! this factor has bitten. Host side only; the caller moves it.
   !
   subroutine basis_subset(this, s0, s1, out)
      class(trc_basis_t), intent(in) :: this
      integer, intent(in) :: s0, s1
      type(trc_basis_t), intent(out) :: out
      integer :: i, n
      n = s1 - s0 + 1
      out%nshell = n
      out%natm = this%natm
      out%maxnp = this%maxnp
      allocate (out%sh_l(n), out%sh_np(n), out%sh_ao(n), out%sh_ao_sph(n))
      allocate (out%sh_e(this%maxnp, n), out%sh_c(this%maxnp, n), out%sh_r(3, n))
      allocate (out%sh_at(n))
      allocate (out%at_z(this%natm), out%at_r(3, this%natm), out%at_elem(this%natm))
      out%sh_l = this%sh_l(s0:s1)
      out%sh_np = this%sh_np(s0:s1)
      out%sh_e = this%sh_e(:, s0:s1)
      out%sh_c = this%sh_c(:, s0:s1)
      out%sh_r = this%sh_r(:, s0:s1)
      out%sh_at = this%sh_at(s0:s1)
      out%at_z = this%at_z
      out%at_r = this%at_r
      out%at_elem = this%at_elem
      ! Same convention as the parent -- an atomic SAD/SADQ subbasis stays
      ! spherical when the molecule it was cut from is, so its own atomic
      ! SCF (trc_sad.F90) runs in the same AO space and the density it
      ! hands back slots straight into the molecular one with no further
      ! conversion at that boundary.
      out%spherical = this%spherical
      out%nao = 0
      out%nao_sph = 0
      do i = 1, n
         out%sh_ao(i) = out%nao + 1
         out%nao = out%nao + ncart(out%sh_l(i))
         out%sh_ao_sph(i) = out%nao_sph + 1
         out%nao_sph = out%nao_sph + 2*out%sh_l(i) + 1
      end do
   end subroutine basis_subset

   subroutine basis_from_libcint(this, atm, natm, bas, nshell, env, at_elem)
      class(trc_basis_t), intent(inout) :: this
      integer,  intent(in) :: natm, nshell
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)
      integer,  intent(in), optional :: at_elem(natm)
      !! Element per atom. libcint's atm carries only the CHARGE (0 on a ghost),
      !! so without this every ghost's element is nint(0) = 0: an element-0
      !! DFT grid and Becke partition at a ghost oxygen. Absent, as before.

      integer :: i, k, ia, np, pe, pc, maxnp
      integer,  allocatable :: l(:), nprim(:)
      real(dp), allocatable :: e(:, :), c(:, :), r(:, :), z(:), rat(:, :)

      maxnp = 1
      do i = 0, nshell - 1
         maxnp = max(maxnp, bas(BAS_SLOTS*i + NPRIM_OF))
      end do

      allocate (l(nshell), nprim(nshell), r(3, nshell))
      allocate (e(maxnp, nshell), c(maxnp, nshell))
      allocate (z(natm), rat(3, natm))
      e = 0.0_dp
      c = 0.0_dp

      do ia = 1, natm
         z(ia) = real(atm(ATM_SLOTS*(ia - 1) + CHARGE_OF), dp)
         rat(:, ia) = env(atm(ATM_SLOTS*(ia - 1) + PTR_COORD): &
                          atm(ATM_SLOTS*(ia - 1) + PTR_COORD) + 2)
      end do

      do i = 1, nshell
         l(i)     = bas(BAS_SLOTS*(i - 1) + ANG_OF)
         nprim(i) = bas(BAS_SLOTS*(i - 1) + NPRIM_OF)
         pe       = bas(BAS_SLOTS*(i - 1) + PTR_EXP)
         pc       = bas(BAS_SLOTS*(i - 1) + PTR_COEFF)
         ia       = bas(BAS_SLOTS*(i - 1) + ATOM_OF)
         r(:, i)  = env(atm(ATM_SLOTS*ia + PTR_COORD): &
                        atm(ATM_SLOTS*ia + PTR_COORD) + 2)
         do k = 1, nprim(i)
            e(k, i) = env(pe + k - 1)
            c(k, i) = env(pc + k - 1)
         end do
      end do

      call this%build(nshell, l, nprim, e, c, r, natm, z, rat, maxnp, at_elem=at_elem)
      deallocate (l, nprim, e, c, r, z, rat)
   end subroutine basis_from_libcint

   !
   ! Parent first, then every allocatable component -- and the reverse on the
   ! way out. A half-mapped derived type does not fail loudly; the mapped
   ! fields work and the forgotten one reads garbage.
   !
   subroutine basis_to_device(this)
      class(trc_basis_t), intent(inout) :: this
      if (this%on_device) return
      !$acc enter data copyin(this)
      !$acc enter data copyin(this%sh_l, this%sh_np, this%sh_ao, this%sh_ao_sph, &
      !$acc                   this%sh_e, this%sh_c, this%sh_r, this%sh_at, &
      !$acc                   this%at_z, this%at_r)
      this%on_device = .true.
   end subroutine basis_to_device

   subroutine basis_release(this)
      class(trc_basis_t), intent(inout) :: this
      if (this%on_device) then
         !$acc exit data delete(this%sh_l, this%sh_np, this%sh_ao, this%sh_ao_sph, &
         !$acc                  this%sh_e, this%sh_c, this%sh_r, this%sh_at, &
         !$acc                  this%at_z, this%at_r)
         !$acc exit data delete(this)
         this%on_device = .false.
      end if
      if (allocated(this%sh_l))  deallocate (this%sh_l)
      if (allocated(this%sh_np)) deallocate (this%sh_np)
      if (allocated(this%sh_ao)) deallocate (this%sh_ao)
      if (allocated(this%sh_ao_sph)) deallocate (this%sh_ao_sph)
      if (allocated(this%sh_e))  deallocate (this%sh_e)
      if (allocated(this%sh_c))  deallocate (this%sh_c)
      if (allocated(this%sh_r))  deallocate (this%sh_r)
      if (allocated(this%sh_at)) deallocate (this%sh_at)
      if (allocated(this%at_z))  deallocate (this%at_z)
      if (allocated(this%at_r))  deallocate (this%at_r)
      if (allocated(this%at_elem)) deallocate (this%at_elem)
      this%nshell = 0
      this%nao = 0
      this%natm = 0
      this%nao_sph = 0
      this%spherical = .false.
   end subroutine basis_release

   subroutine pairs_build(this, b, thresh)
      class(trc_pairlist_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: thresh

      integer  :: i, j, ki, kj, n, pass, np, d
      real(dp) :: al, bet, zeta, xi, r2, kab, sb, vb, zsum, dx, dy, dz

      call this%release()
      this%nshell = b%nshell
      this%thresh = thresh

      zsum = 0.0_dp
      do i = 1, b%natm
         zsum = zsum + abs(b%at_z(i))
      end do
      if (zsum <= 0.0_dp) zsum = 1.0_dp

      ! Two passes: count, then fill. The pair list is O(N^2) in shells and
      ! the bound is cheap, so counting twice beats guessing a capacity.
      do pass = 1, 2
         n = 0
         do i = 1, b%nshell
            do j = 1, i
               dx = b%sh_r(1, i) - b%sh_r(1, j)
               dy = b%sh_r(2, i) - b%sh_r(2, j)
               dz = b%sh_r(3, i) - b%sh_r(3, j)
               r2 = dx*dx + dy*dy + dz*dz
               sb = 0.0_dp
               vb = 0.0_dp
               do ki = 1, b%sh_np(i)
                  al = b%sh_e(ki, i)
                  do kj = 1, b%sh_np(j)
                     bet = b%sh_e(kj, j)
                     zeta = al + bet
                     xi = al*bet/zeta
                     kab = abs(b%sh_c(ki, i)*b%sh_c(kj, j))*exp(-xi*r2)
                     sb = sb + kab*(PI/zeta)**1.5_dp
                     vb = vb + kab*(2.0_dp*PI/zeta)
                  end do
               end do
               if (max(sb, vb*zsum) <= thresh) cycle
               n = n + 1
               if (pass == 2) then
                  this%p_i(n) = i
                  this%p_j(n) = j
                  this%p_s(n) = sb
                  this%p_v(n) = vb
                  this%p_n(n) = b%sh_np(i)*b%sh_np(j)
               end if
            end do
         end do
         if (pass == 1) then
            this%npair = n
            allocate (this%p_i(n), this%p_j(n), this%p_s(n), this%p_v(n))
            allocate (this%p_n(n))
         end if
      end do

      ! second stage: the primitive pairs themselves, prefix-summed
      allocate (this%pp_off(this%npair), this%pp_n(this%npair))
      this%npp = 0
      do n = 1, this%npair
         this%pp_off(n) = this%npp
         this%pp_n(n) = this%p_n(n)
         this%npp = this%npp + this%p_n(n)
      end do
      allocate (this%pp_z(this%npp), this%pp_k(this%npp))
      allocate (this%pp_p(3, this%npp), this%pp_a(3, this%npp))

      np = 0
      do n = 1, this%npair
         i = this%p_i(n)
         j = this%p_j(n)
         r2 = 0.0_dp
         do d = 1, 3
            r2 = r2 + (b%sh_r(d, i) - b%sh_r(d, j))**2
         end do
         do ki = 1, b%sh_np(i)
            al = b%sh_e(ki, i)
            do kj = 1, b%sh_np(j)
               bet = b%sh_e(kj, j)
               zeta = al + bet
               xi = al*bet/zeta
               np = np + 1
               this%pp_z(np) = zeta
               this%pp_k(np) = b%sh_c(ki, i)*b%sh_c(kj, j)*exp(-xi*r2)
               do d = 1, 3
                  this%pp_p(d, np) = (al*b%sh_r(d, i) + bet*b%sh_r(d, j))/zeta
                  this%pp_a(d, np) = this%pp_p(d, np) - b%sh_r(d, i)
               end do
            end do
         end do
      end do
   end subroutine pairs_build

   subroutine pairs_to_device(this)
      class(trc_pairlist_t), intent(inout) :: this
      if (this%on_device) return
      !$acc enter data copyin(this)
      !$acc enter data copyin(this%p_i, this%p_j, this%p_s, this%p_v, &
      !$acc                   this%pp_off, this%pp_n, this%pp_z, this%pp_k, &
      !$acc                   this%pp_p, this%pp_a)
      this%on_device = .true.
   end subroutine pairs_to_device

   subroutine pairs_release(this)
      class(trc_pairlist_t), intent(inout) :: this
      if (this%on_device) then
         !$acc exit data delete(this%p_i, this%p_j, this%p_s, this%p_v, &
         !$acc                  this%pp_off, this%pp_n, this%pp_z, this%pp_k, &
         !$acc                  this%pp_p, this%pp_a)
         !$acc exit data delete(this)
         this%on_device = .false.
      end if
      if (allocated(this%p_i)) deallocate (this%p_i)
      if (allocated(this%p_j)) deallocate (this%p_j)
      if (allocated(this%p_s)) deallocate (this%p_s)
      if (allocated(this%p_v)) deallocate (this%p_v)
      if (allocated(this%p_n)) deallocate (this%p_n)
      if (allocated(this%pp_off)) deallocate (this%pp_off)
      if (allocated(this%pp_n)) deallocate (this%pp_n)
      if (allocated(this%pp_z)) deallocate (this%pp_z)
      if (allocated(this%pp_k)) deallocate (this%pp_k)
      if (allocated(this%pp_p)) deallocate (this%pp_p)
      if (allocated(this%pp_a)) deallocate (this%pp_a)
      this%npair = 0
      this%npp = 0
   end subroutine pairs_release

end module trc_basis_types
