!
! The consumer-facing API: a device-resident basis object and batched drivers.
!
! WHAT THIS IS FOR
! ----------------
! Everything below the kernels is per-class `!$acc routine seq` code that
! computes ONE shell block. That is the right shape for a generator and the
! wrong shape for a caller. This module is the boundary: a caller hands over a
! basis once, the data goes to the device once and stays there, and each
! driver is a single `do concurrent` over shell pairs.
!
! SHAPE OF THE INTERFACE
! ----------------------
! The basis is accepted in libcint's `atm`/`bas`/`env` layout, because that is
! what mqc already has -- it carries a libcint backend, so no caller-side
! conversion is needed to try this. It is converted ONCE into a
! structure-of-arrays container that the kernels can read coalesced, rather
! than being walked in that packed integer form on the device.
!
!     type(trc_basis_t) :: b
!     call b%from_libcint(atm, natm, bas, nshell, env)
!     call b%to_device()
!     call trc_1e(b, smat, tmat, vmat)      ! nao x nao each
!     call b%release()
!
! WHY ONE CALL FOR S, T AND V
! ---------------------------
! They share every shell-pair quantity -- zeta, P, the Gaussian product
! prefactor -- and a caller forming a core Hamiltonian wants all three. Asking
! for them separately would read the primitive data three times to save
! nothing.
!
! NO ATOMICS HERE
! ---------------
! Distinct shell pairs write disjoint blocks of the output matrices, so each
! thread owns its output outright. That is the one respect in which the
! one-electron drivers are simpler than the Fock build, which has to scatter
! into six overlapping blocks per quartet.
!
module trc_api
   use trc_boys, only: dp, boys_init
   use trc_tables, only: LMAX
   use trc_1e_kernels, only: one_e_dispatch, ONE_E_LMAX
   use trc_mult_kernels, only: multipole_dispatch, TRC_NMULT, MULT_LMAX
   use trc_df_kernels, only: df_2c_dispatch, df_3c_dispatch, &
                               DF_LMAX, DF_LMAX_AUX
   implicit none
   private
   public :: trc_basis_t, trc_pairlist_t
   public :: trc_1e, trc_multipoles, trc_df_2c, trc_df_3c
   public :: TRC_NMULT

   ! libcint's packed layout, so a caller can pass its arrays straight in.
   integer, parameter :: ATM_SLOTS = 6, BAS_SLOTS = 8
   integer, parameter :: CHARGE_OF = 0, PTR_COORD = 1
   integer, parameter :: ATOM_OF = 0, ANG_OF = 1, NPRIM_OF = 2
   integer, parameter :: PTR_EXP = 5, PTR_COEFF = 6

   !
   ! Angular-momentum limits, taken FROM THE GENERATORS rather than restated.
   !
   ! The four-centre ERIs and everything else have different ceilings on
   ! purpose. A four-centre class carries momentum on all four centres, so its
   ! kernels grow as L^4 and (dd|dd) is already 20k lines; the one-electron and
   ! density-fitting classes carry it on two or three, so f and g are
   ! affordable there while they would not be for the ERIs. Splitting the two
   ! is what lets the fitting basis outrank the orbital basis, which is how
   ! fitting bases are normally built.
   !
   ! TRC_MAXL is the limit for THIS module's paths (1e and DF). The
   ! four-centre limit is TRC_LMAX, a preprocessor symbol consumed by
   ! trc_binkernel and its generated kernels; the two do not have to agree.
   !
   integer, parameter, public :: TRC_MAXL = &
      max(ONE_E_LMAX, max(DF_LMAX, DF_LMAX_AUX))

   !
   ! COMPILE-TIME ASSERTION: the one-electron kernels must reach as high as the
   ! four-centre ones.
   !
   ! "The two do not have to agree" above is true only in the direction where
   ! 1e is HIGHER. The other direction is not a limitation, it is a wrong
   ! answer: a build whose four-centre path outranks ONE_E_LMAX accepts a basis
   ! it cannot do one-electron integrals for, and builds H from whatever
   ! one_e_dispatch returns for a key it has no case for. That is a converged
   ! SCF on a wrong Hamiltonian, with no diagnostic anywhere.
   !
   ! It holds today because both ceilings are 2. It is asserted so that raising
   ! one without the other stops the build instead of the science.
   !
   ! Against trc_tables' LMAX and not the TRC_LMAX preprocessor symbol: CMake
   ! defines TRC_LMAX and fpm does not, so referring to it here compiles under
   ! one build system and not the other. trc_tables already carries the
   ! `#ifndef TRC_LMAX / #define TRC_LMAX 2` default and exports the result as
   ! a parameter, which is the value both builds actually use.
   !
   ! Division rather than a zero-size array because a negative extent is legal
   ! Fortran and silently gives an empty array, whereas 1/0 in a constant
   ! expression is a hard error in gfortran, ifx and nvfortran alike.
   !
   integer, parameter, private :: ASSERT_1E_COVERS_ERI = &
      1/merge(1, 0, ONE_E_LMAX >= LMAX)
   !! Largest Cartesian block a shell pair can produce here. Fixed at compile
   !! time because `do concurrent` locals must be, and because an in-kernel
   !! allocate would stall the device.
   integer, parameter :: NC_MAX = (TRC_MAXL + 1)*(TRC_MAXL + 2)/2
   integer, parameter :: NBLK_MAX = NC_MAX*NC_MAX
   integer, parameter :: NBLK3_MAX = NC_MAX*NC_MAX*NC_MAX
   real(dp), parameter :: PI = 3.14159265358979323846_dp

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
      real(dp), allocatable :: at_z(:)      !! (natm) charges
      real(dp), allocatable :: at_r(:, :)   !! (3, natm) positions
      logical :: on_device = .false.
   contains
      !! The primary constructor: flat arrays, one entry per shell. A caller
      !! that has a basis set has these; only a caller that already links
      !! libcint has `atm`/`bas`/`env`.
      procedure :: build        => basis_build
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

   pure integer function ncart(l)
      integer, intent(in) :: l
      ncart = (l + 1)*(l + 2)/2
   end function ncart

   !
   ! libcint applies this per shell OUTSIDE `env`, so a container built from
   ! `env` alone is wrong by it. Folding it into sh_c here means the kernels
   ! never have to know.
   !
   pure real(dp) function common_fac_sp(l)
      integer, intent(in) :: l
      select case (l)
      case (0); common_fac_sp = 0.282094791773878143_dp
      case (1); common_fac_sp = 0.488602511902919921_dp
      case default; common_fac_sp = 1.0_dp
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
                          natm, at_z, at_r, maxnp)
      class(trc_basis_t), intent(inout) :: this
      integer,  intent(in) :: nshell, natm, maxnp
      integer,  intent(in) :: sh_l(nshell), sh_np(nshell)
      real(dp), intent(in) :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in) :: sh_r(3, nshell)
      real(dp), intent(in) :: at_z(natm), at_r(3, natm)

      integer :: i, k

      call this%release()
      call boys_init()

      this%nshell = nshell
      this%natm   = natm
      this%maxnp  = maxnp

      allocate (this%sh_l(nshell), this%sh_np(nshell), this%sh_ao(nshell))
      allocate (this%sh_e(maxnp, nshell), this%sh_c(maxnp, nshell))
      allocate (this%sh_r(3, nshell))
      allocate (this%at_z(natm), this%at_r(3, natm))
      this%sh_e = 0.0_dp; this%sh_c = 0.0_dp

      this%at_z = at_z
      this%at_r = at_r

      this%nao = 0
      do i = 1, nshell
         this%sh_l(i)  = sh_l(i)
         this%sh_np(i) = sh_np(i)
         this%sh_r(:, i) = sh_r(:, i)
         do k = 1, sh_np(i)
            this%sh_e(k, i) = sh_e(k, i)
            this%sh_c(k, i) = sh_c(k, i)*common_fac_sp(sh_l(i))
         end do
         this%sh_ao(i) = this%nao + 1
         this%nao = this%nao + ncart(sh_l(i))
      end do
   end subroutine basis_build

   !
   ! Convenience adapter for callers that already hold libcint's layout.
   ! NOT the primary interface -- it unpacks into raw arrays and delegates to
   ! `build`, so there is one implementation of the import rules (notably the
   ! common_fac_sp folding) rather than two that can drift.
   !
   subroutine basis_from_libcint(this, atm, natm, bas, nshell, env)
      class(trc_basis_t), intent(inout) :: this
      integer,  intent(in) :: natm, nshell
      integer,  intent(in) :: atm(0:), bas(0:)
      real(dp), intent(in) :: env(0:)

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
      e = 0.0_dp; c = 0.0_dp

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

      call this%build(nshell, l, nprim, e, c, r, natm, z, rat, maxnp)
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
      !$acc enter data copyin(this%sh_l, this%sh_np, this%sh_ao, &
      !$acc                   this%sh_e, this%sh_c, this%sh_r, &
      !$acc                   this%at_z, this%at_r)
      this%on_device = .true.
   end subroutine basis_to_device

   subroutine basis_release(this)
      class(trc_basis_t), intent(inout) :: this
      if (this%on_device) then
         !$acc exit data delete(this%sh_l, this%sh_np, this%sh_ao, &
         !$acc                  this%sh_e, this%sh_c, this%sh_r, &
         !$acc                  this%at_z, this%at_r)
         !$acc exit data delete(this)
         this%on_device = .false.
      end if
      if (allocated(this%sh_l))  deallocate (this%sh_l)
      if (allocated(this%sh_np)) deallocate (this%sh_np)
      if (allocated(this%sh_ao)) deallocate (this%sh_ao)
      if (allocated(this%sh_e))  deallocate (this%sh_e)
      if (allocated(this%sh_c))  deallocate (this%sh_c)
      if (allocated(this%sh_r))  deallocate (this%sh_r)
      if (allocated(this%at_z))  deallocate (this%at_z)
      if (allocated(this%at_r))  deallocate (this%at_r)
      this%nshell = 0; this%nao = 0; this%natm = 0
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
               sb = 0.0_dp; vb = 0.0_dp
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
                  this%p_i(n) = i; this%p_j(n) = j
                  this%p_s(n) = sb; this%p_v(n) = vb
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
         i = this%p_i(n); j = this%p_j(n)
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
      this%npair = 0; this%npp = 0
   end subroutine pairs_release

   !
   ! Cartesian multipole moments about `orig`, through the octupole.
   !
   ! (nao, nao, 39): 3 dipole, 9 quadrupole, 27 octupole, full tensors in
   ! libcint's order with the last index fastest -- what makeEFP's distributed
   ! multipoles consume.
   !
   ! A separate call from `trc_1e` even though the two share 1-D overlap
   ! tables. Fusing them was tried and put 39 components per class into one
   ! device routine, which NVVM chewed on for over thirty minutes at LMAX=2
   ! against two and a half apart. Both are computed once per geometry, so the
   ! shared tables were never worth it.
   !
   subroutine trc_multipoles(b, p, orig, mmat)
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: p
      real(dp), intent(in)  :: orig(3)
      real(dp), intent(out) :: mmat(b%nao, b%nao, TRC_NMULT)
      call mult_driver(b%nshell, b%nao, b%maxnp, b%sh_l, b%sh_np, b%sh_ao, &
                       b%sh_e, b%sh_c, b%sh_r, p%npair, p%p_i, p%p_j, orig, mmat)
   end subroutine trc_multipoles

   subroutine mult_driver(nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                          sh_e, sh_c, sh_r, npair, p_i, p_j, orig, mmat)
      integer,  intent(in)  :: nshell, nao, maxnp, npair
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell), orig(3)
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      real(dp), intent(out) :: mmat(nao, nao, TRC_NMULT)

      integer :: ij, i, j, k

      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      do concurrent(k=1:TRC_NMULT, j=1:nao, i=1:nao)
         mmat(i, j, k) = 0.0_dp
      end do

      do concurrent(ij = 1:npair)
         call mult_item(p_i(ij), p_j(ij), nshell, nao, maxnp, sh_l, sh_np, &
                        sh_ao, sh_e, sh_c, sh_r, orig, mmat)
      end do
   end subroutine mult_driver

   pure subroutine mult_item(ii, jj, nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                             sh_e, sh_c, sh_r, orig, mmat)
      !$acc routine seq
      integer,  intent(in)    :: ii, jj, nshell, nao, maxnp
      integer,  intent(in)    :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)    :: sh_r(3, nshell), orig(3)
      real(dp), intent(inout) :: mmat(nao, nao, TRC_NMULT)

      integer  :: li, lj, ni, nj, a, c, mu, nu, k, ic
      !
      ! ONE DIMENSION, INDEXED BY HAND. The generated kernel declares its
      ! output `mout(ni*nj, 39)` with the ACTUAL block size of its class,
      ! and reaches this routine through a `mout(*)` dummy -- so it writes
      ! component ic starting at (ic-1)*ni*nj. Declaring the buffer here as
      ! `mblk(NBLK_MAX, TRC_NMULT)` reads it back with stride NBLK_MAX
      ! instead, which agrees only for the one class whose block happens to
      ! be NBLK_MAX and silently scrambles every other one. That bug cost a
      ! day: the dipole x component is the first ni*nj elements and so came
      ! out right, which made it look like an ordering problem in the
      ! generator rather than a stride mismatch here.
      !
      real(dp) :: mblk(NBLK_MAX*TRC_NMULT)

      li = sh_l(ii); lj = sh_l(jj)
      ni = (li + 1)*(li + 2)/2
      nj = (lj + 1)*(lj + 2)/2
      call multipole_dispatch(li, lj, sh_np(ii), sh_np(jj), &
                              sh_e(:, ii), sh_c(:, ii), &
                              sh_e(:, jj), sh_c(:, jj), &
                              sh_r(:, ii), sh_r(:, jj), orig, mblk)
      k = 0
      do c = 0, nj - 1
         nu = sh_ao(jj) + c
         do a = 0, ni - 1
            mu = sh_ao(ii) + a
            k = k + 1
            do ic = 1, TRC_NMULT
               ! the multipole matrices are symmetric: the operator is
               ! multiplicative and carries no derivative
               mmat(mu, nu, ic) = mblk(k + (ic - 1)*ni*nj)
               mmat(nu, mu, ic) = mblk(k + (ic - 1)*ni*nj)
            end do
         end do
      end do
   end subroutine mult_item

   !> Overlap, kinetic and nuclear attraction, all three at once.
   subroutine trc_1e(b, p, smat, tmat, vmat)
      type(trc_basis_t), intent(in)  :: b
      type(trc_pairlist_t), intent(in) :: p
      real(dp), intent(out) :: smat(b%nao, b%nao)
      real(dp), intent(out) :: tmat(b%nao, b%nao)
      real(dp), intent(out) :: vmat(b%nao, b%nao)
      call one_e_driver(b%nshell, b%nao, b%natm, b%maxnp, b%sh_l, b%sh_np, &
                        b%sh_ao, b%sh_e, b%sh_c, b%sh_r, b%at_z, b%at_r, &
                        p%npair, p%p_i, p%p_j, smat, tmat, vmat)
   end subroutine trc_1e

   !
   ! Explicit-shape dummies with the extents passed in, and the extents
   ! declared first. An assumed-shape dummy makes nvfortran walk a descriptor
   ! per launch; reaching through b%component inside the loop does the same.
   !
   subroutine one_e_driver(nshell, nao, natm, maxnp, sh_l, sh_np, sh_ao, &
                           sh_e, sh_c, sh_r, at_z, at_r, &
                           npair, p_i, p_j, smat, tmat, vmat)
      integer,  intent(in)  :: nshell, nao, natm, maxnp, npair
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell)
      real(dp), intent(in)  :: at_z(natm), at_r(3, natm)
      real(dp), intent(out) :: smat(nao, nao), tmat(nao, nao), vmat(nao, nao)

      integer :: ij, i, j

      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      do concurrent(j=1:nao, i=1:nao)
         smat(i, j) = 0.0_dp
         tmat(i, j) = 0.0_dp
         vmat(i, j) = 0.0_dp
      end do

      ! One thread per SURVIVING shell pair. The list is lower-triangular and
      ! these matrices are symmetric, so each thread writes its block and the
      ! transpose -- which is why the earlier full-square version existed. The
      ! pair list changes the arithmetic: with screening, walking the square
      ! costs a factor of two on top of computing pairs that are known to be
      ! negligible.
      !
      ! The body lives in a `pure` item routine rather than inline, for two
      ! reasons: nvfortran rejects a BLOCK construct inside `do concurrent`,
      ! and declaring the scratch alongside the loop would make it SHARED
      ! across threads -- the failure that silently corrupted the four-centre
      ! per-class kernels once already.
      do concurrent(ij = 1:npair)
         call one_e_item(p_i(ij), p_j(ij), nshell, nao, natm, maxnp, sh_l, &
                         sh_np, sh_ao, sh_e, sh_c, sh_r, at_z, at_r, &
                         smat, tmat, vmat)
      end do
   end subroutine one_e_driver

   pure subroutine one_e_item(ii, jj, nshell, nao, natm, maxnp, sh_l, sh_np, &
                              sh_ao, sh_e, sh_c, sh_r, at_z, at_r, &
                              smat, tmat, vmat)
      !$acc routine seq
      integer,  intent(in)    :: ii, jj, nshell, nao, natm, maxnp
      integer,  intent(in)    :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)    :: sh_r(3, nshell)
      real(dp), intent(in)    :: at_z(natm), at_r(3, natm)
      real(dp), intent(inout) :: smat(nao, nao), tmat(nao, nao), vmat(nao, nao)

      integer  :: li, lj, ni, nj, a, c, mu, nu, k
      real(dp) :: sblk(NBLK_MAX), tblk(NBLK_MAX), vblk(NBLK_MAX)

      li = sh_l(ii); lj = sh_l(jj)
      ni = (li + 1)*(li + 2)/2
      nj = (lj + 1)*(lj + 2)/2
      call one_e_dispatch(li, lj, sh_np(ii), sh_np(jj), &
                          sh_e(:, ii), sh_c(:, ii), &
                          sh_e(:, jj), sh_c(:, jj), &
                          sh_r(:, ii), sh_r(:, jj), &
                          natm, at_z, at_r, sblk, tblk, vblk)
      k = 0
      do c = 0, nj - 1
         nu = sh_ao(jj) + c
         do a = 0, ni - 1
            mu = sh_ao(ii) + a
            k = k + 1
            smat(mu, nu) = sblk(k)
            tmat(mu, nu) = tblk(k)
            vmat(mu, nu) = vblk(k)
            ! symmetric counterpart; on the diagonal block (ii == jj) this
            ! rewrites the same element with the same value, which is why it
            ! needs no guard
            smat(nu, mu) = sblk(k)
            tmat(nu, mu) = tblk(k)
            vmat(nu, mu) = vblk(k)
         end do
      end do
   end subroutine one_e_item

   !> (P|Q) over an auxiliary basis.
   subroutine trc_df_2c(aux, jmat)
      type(trc_basis_t), intent(in) :: aux
      real(dp), intent(out) :: jmat(aux%nao, aux%nao)
      call df_2c_driver(aux%nshell, aux%nao, aux%maxnp, aux%sh_l, aux%sh_np, &
                        aux%sh_ao, aux%sh_e, aux%sh_c, aux%sh_r, jmat)
   end subroutine trc_df_2c

   subroutine df_2c_driver(nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                           sh_e, sh_c, sh_r, jmat)
      integer,  intent(in)  :: nshell, nao, maxnp
      integer,  intent(in)  :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)  :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)  :: sh_r(3, nshell)
      real(dp), intent(out) :: jmat(nao, nao)

      integer :: ij, i, j

      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      do concurrent(j=1:nao, i=1:nao)
         jmat(i, j) = 0.0_dp
      end do

      do concurrent(ij = 0:nshell*nshell - 1)
         call df_2c_item(ij, nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                         sh_e, sh_c, sh_r, jmat)
      end do
   end subroutine df_2c_driver

   pure subroutine df_2c_item(ij, nshell, nao, maxnp, sh_l, sh_np, sh_ao, &
                              sh_e, sh_c, sh_r, jmat)
      !$acc routine seq
      integer,  intent(in)    :: ij, nshell, nao, maxnp
      integer,  intent(in)    :: sh_l(nshell), sh_np(nshell), sh_ao(nshell)
      real(dp), intent(in)    :: sh_e(maxnp, nshell), sh_c(maxnp, nshell)
      real(dp), intent(in)    :: sh_r(3, nshell)
      real(dp), intent(inout) :: jmat(nao, nao)

      integer  :: ii, jj, li, lj, ni, nj, a, c, mu, nu, k
      real(dp) :: blk(NBLK_MAX)

      ii = ij/nshell + 1
      jj = mod(ij, nshell) + 1
      li = sh_l(ii); lj = sh_l(jj)
      ni = (li + 1)*(li + 2)/2
      nj = (lj + 1)*(lj + 2)/2
      call df_2c_dispatch(li, lj, sh_np(ii), sh_e(:, ii), sh_c(:, ii), &
                          sh_r(:, ii), sh_np(jj), sh_e(:, jj), &
                          sh_c(:, jj), sh_r(:, jj), blk)
      k = 0
      do c = 0, nj - 1
         nu = sh_ao(jj) + c
         do a = 0, ni - 1
            mu = sh_ao(ii) + a
            k = k + 1
            jmat(mu, nu) = blk(k)
         end do
      end do
   end subroutine df_2c_item

   !> (mu nu|P), the three-index tensor density fitting is built on.
   subroutine trc_df_3c(b, p, aux, tens)
      type(trc_basis_t), intent(in) :: b, aux
      type(trc_pairlist_t), intent(in) :: p
      real(dp), intent(out) :: tens(b%nao, b%nao, aux%nao)
      call df_3c_driver(b%nshell, b%nao, b%maxnp, b%sh_l, b%sh_ao, b%sh_r, &
                        p%npair, p%npp, p%p_i, p%p_j, p%pp_off, p%pp_n, &
                        p%pp_z, p%pp_k, p%pp_p, p%pp_a, &
                        aux%nshell, aux%nao, aux%maxnp, aux%sh_l, aux%sh_np, &
                        aux%sh_ao, aux%sh_e, aux%sh_c, aux%sh_r, tens)
   end subroutine trc_df_3c

   subroutine df_3c_driver(nsh, nao, mnp, sh_l, sh_ao, sh_r, &
                           npair, npp, p_i, p_j, pp_off, pp_n, &
                           pp_z, pp_k, pp_p, pp_a, &
                           nsx, nax, mnx, ax_l, ax_np, ax_ao, ax_e, ax_c, ax_r, &
                           tens)
      integer,  intent(in)  :: nsh, nao, mnp, npair, npp, nsx, nax, mnx
      integer,  intent(in)  :: sh_l(nsh), sh_ao(nsh)
      real(dp), intent(in)  :: sh_r(3, nsh)
      integer,  intent(in)  :: p_i(npair), p_j(npair)
      integer,  intent(in)  :: pp_off(npair), pp_n(npair)
      real(dp), intent(in)  :: pp_z(npp), pp_k(npp), pp_p(3, npp), pp_a(3, npp)
      integer,  intent(in)  :: ax_l(nsx), ax_np(nsx), ax_ao(nsx)
      real(dp), intent(in)  :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      real(dp), intent(out) :: tens(nao, nao, nax)

      integer(kind=8) :: t
      integer :: i, j, k

      ! `do concurrent` rather than an OpenACC parallel loop. The arrays are
      ! already device-resident, so the standard construct is enough and the
      ! directive was only ever spelling out what it already says.
      do concurrent(k=1:nax, j=1:nao, i=1:nao)
         tens(i, j, k) = 0.0_dp
      end do

      ! One thread per (SURVIVING shell pair, auxiliary shell). The pair index
      ! runs fastest so that neighbouring threads share the auxiliary shell,
      ! and therefore its primitive data, across a warp.
      do concurrent(t = 0:int(npair, 8)*nsx - 1)
         call df_3c_item(t, nsh, nao, sh_l, sh_ao, sh_r, npair, npp, &
                         p_i, p_j, pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, &
                         nsx, nax, mnx, ax_l, ax_np, ax_ao, ax_e, ax_c, ax_r, &
                         tens)
      end do
   end subroutine df_3c_driver

   pure subroutine df_3c_item(t, nsh, nao, sh_l, sh_ao, sh_r, npair, npp, &
                              p_i, p_j, pp_off, pp_n, pp_z, pp_k, pp_p, pp_a, &
                              nsx, nax, mnx, ax_l, ax_np, ax_ao, ax_e, ax_c, &
                              ax_r, tens)
      !$acc routine seq
      integer(kind=8), intent(in) :: t
      integer,  intent(in)    :: nsh, nao, npair, npp, nsx, nax, mnx
      integer,  intent(in)    :: sh_l(nsh), sh_ao(nsh)
      real(dp), intent(in)    :: sh_r(3, nsh)
      integer,  intent(in)    :: p_i(npair), p_j(npair)
      integer,  intent(in)    :: pp_off(npair), pp_n(npair)
      real(dp), intent(in)    :: pp_z(npp), pp_k(npp)
      real(dp), intent(in)    :: pp_p(3, npp), pp_a(3, npp)
      integer,  intent(in)    :: ax_l(nsx), ax_np(nsx), ax_ao(nsx)
      real(dp), intent(in)    :: ax_e(mnx, nsx), ax_c(mnx, nsx), ax_r(3, nsx)
      real(dp), intent(inout) :: tens(nao, nao, nax)

      integer :: ip, kk, ii, jj, li, lj, lk, ni, nj, nk
      integer :: a, c, e, mu, nu, pq, m, o
      real(dp) :: blk(NBLK3_MAX), abx, aby, abz

      ip = int(mod(t, int(npair, 8))) + 1
      kk = int(t/npair) + 1
      ii = p_i(ip); jj = p_j(ip)
      li = sh_l(ii); lj = sh_l(jj); lk = ax_l(kk)
      ni = (li + 1)*(li + 2)/2
      nj = (lj + 1)*(lj + 2)/2
      nk = (lk + 1)*(lk + 2)/2
      abx = sh_r(1, ii) - sh_r(1, jj)
      aby = sh_r(2, ii) - sh_r(2, jj)
      abz = sh_r(3, ii) - sh_r(3, jj)
      o = pp_off(ip)

      call df_3c_dispatch(li, lj, lk, pp_n(ip), &
                          pp_z(o + 1:o + pp_n(ip)), pp_k(o + 1:o + pp_n(ip)), &
                          pp_p(:, o + 1:o + pp_n(ip)), &
                          pp_a(:, o + 1:o + pp_n(ip)), &
                          abx, aby, abz, &
                          ax_np(kk), ax_e(:, kk), ax_c(:, kk), ax_r(:, kk), blk)

      ! The pair list is lower-triangular in (ii, jj), and (mu nu|P) is
      ! symmetric in mu <-> nu, so each thread fills both. On a diagonal
      ! shell pair the second write repeats the first with the same value.
      m = 0
      do e = 0, nk - 1
         pq = ax_ao(kk) + e
         do c = 0, nj - 1
            nu = sh_ao(jj) + c
            do a = 0, ni - 1
               mu = sh_ao(ii) + a
               m = m + 1
               tens(mu, nu, pq) = blk(m)
               tens(nu, mu, pq) = blk(m)
            end do
         end do
      end do
   end subroutine df_3c_item

end module trc_api
