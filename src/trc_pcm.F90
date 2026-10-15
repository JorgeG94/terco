!
! A polarisable continuum: the cavity, the surface charges, and the
! one-electron operator they add to the Fock matrix.
!
! The smooth switching/Gaussian (SWIG) discretisation of Lange and Herbert
! [J. Chem. Phys. 133, 244111 (2010)], solved as C-PCM or IEF-PCM, following
! pyscf's `pyscf.solvent.pcm` term for term: the same Lebedev points per
! sphere, the same switching function, per-point Gaussian exponents from
! the fitted xi of Table II in [J. Chem. Phys. 122, 194110 (2005)], and the
! same S and D matrices. Ported from metalquicha's mqc_czt_pcm, which
! follows pyscf the same way; what changed on the way is where the
! solute-surface integrals come from (terco's own kernels, trc_pcm_potential
! and trc_pcm_fock, rather than libcint's int1e_grids) and that the matrix
! builds are `do concurrent`.
!
! WHAT A SURFACE POINT IS
! -----------------------
! Not a point charge but a normalised spherical Gaussian of exponent
! zeta**2. That keeps every matrix element finite as two points approach,
! and the cavity smooth as a point fades through the switching region. Its
! potential at distance d is erf(zeta d)/d, which is what the integral
! kernels evaluate against the basis and what v_nuc evaluates against the
! nuclei.
!
! THE MODELS
! ----------
!     C-PCM:   K = S                          R = -f_eps I
!     IEF-PCM: K = S - f_eps/(2 pi) D A S     R = -f_eps (I - D A/(2 pi))
!
! with f_eps = (eps - 1)/eps for C-PCM and (eps - 1)/(eps + 1) for IEF-PCM,
! A the point areas, S carrying erf(zeta_ij r)/r off the diagonal and D its
! normal derivative (eqs. 3.35-3.43 of the SWIG paper as pyscf writes them).
! The charges solve K q = R v for the solute potential v on the surface;
! for IEF-PCM's nonsymmetric K the transposed system is solved too and the
! two averaged, as pyscf does, which is what makes E = q.v/2 a variational
! functional of the density. The Fock matrix takes the FULL potential of
! the charges (the half does not survive the derivative).
!
! THE SOLVE, AND WHAT IT COSTS
! ----------------------------
! Dense: K is formed once per geometry and LU-factored, and each SCF
! iteration is then two triangular solves. All of it through trc_linalg --
! cuSOLVER and cuBLAS on matrices that stay on the device under the GPU
! build, LAPACK on the host otherwise -- so an iteration moves n numbers
! across the bus, not n^2, and never meets the host BLAS's threading
! (threaded MKL and NVHPC's OpenMP runtime do not share a process). That
! is O(n^2) memory
! in the number of surface points n -- 8 n^2 bytes for C-PCM's K, three
! such matrices while IEF-PCM's are being formed -- which is fine to a few
! thousand points per molecule (a 100-atom solute at 302 points per sphere
! is ~10^4 points, ~1 GB) and is the wall for large clusters. Everything
! that touches K goes through `solve` and `build_operators`, so a
! matrix-free K q (recomputing S's elements, as the builds below already
! do one at a time) under a conjugate-gradient solve can replace both
! without the SCF noticing.
!
module trc_pcm
   use trc_boys, only: dp
   use trc_linalg, only: trc_linalg_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_lebedev, only: lebedev_grid, lebedev_is_available
   use trc_pcm_radii, only: pcm_cavity_radius
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_pcm_potential, trc_pcm_fock, PCM_LMAX
   use pic_mpi_lib, only: comm_t, allreduce, bcast, MPI_SUM
   implicit none
   private

   public :: trc_pcm_t, pcm_surface_t
   public :: build_pcm_surface, pcm_sd_matrices, pcm_switch_h, pcm_model_from_name
   public :: PCM_MODEL_NONE, PCM_MODEL_CPCM, PCM_MODEL_IEFPCM
   public :: PCM_DEFAULT_LEBEDEV, PCM_XI_ORDERS

   integer, parameter :: PCM_MODEL_NONE = 0
   integer, parameter :: PCM_MODEL_CPCM = 1
      !! conductor-like: K = S, scaled by (eps-1)/eps
   integer, parameter :: PCM_MODEL_IEFPCM = 2
      !! integral equation formalism: the D terms, scaled by (eps-1)/(eps+1)
   integer, parameter :: PCM_DEFAULT_LEBEDEV = 302
      !! points per sphere; pyscf's lebedev_order = 29

   real(dp), parameter :: PI = 3.14159265358979323846_dp

   ! The Lebedev orders the SWIG exponents are fitted for, and the fitted
   ! values: Table II of [J. Chem. Phys. 122, 194110 (2005)], as pyscf's
   ! `pcm.XI` carries them. terco's Lebedev tables also have 74, 230 and
   ! 266, which the paper does not tabulate; those are refused rather than
   ! interpolated, since a wrong exponent smooths the cavity by the wrong
   ! amount and changes the energy without failing.
   integer, parameter :: N_XI = 29
   integer, parameter :: PCM_XI_ORDERS(N_XI) = [ &
                         6, 14, 26, 38, 50, 86, 110, 146, 170, 194, 302, 350, 434, 590, &
                         770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, &
                         4334, 4802, 5294, 5810]
   real(dp), parameter :: XI_VALUES(N_XI) = [ &
                          4.84566077868_dp, 4.86458714334_dp, 4.85478226219_dp, 4.90105812685_dp, &
                          4.89250673295_dp, 4.89741372580_dp, 4.90101060987_dp, 4.89825187392_dp, &
                          4.90685517725_dp, 4.90337644248_dp, 4.90498088169_dp, 4.86879474832_dp, &
                          4.90567349080_dp, 4.90624071359_dp, 4.90656435779_dp, 4.90685167998_dp, &
                          4.90704098216_dp, 4.90721023869_dp, 4.90733270691_dp, 4.90744499142_dp, &
                          4.90753082825_dp, 4.90760972766_dp, 4.90767282394_dp, 4.90773141371_dp, &
                          4.90777965981_dp, 4.90782469526_dp, 4.90749125553_dp, 4.90762073452_dp, &
                          4.90792902522_dp]

   real(dp), parameter :: BURIED_CUTOFF = 1.0e-16_dp
      !! a point whose weight times switching value is at or below this is
      !! buried in a neighbouring sphere and dropped, as pyscf drops it
   real(dp), parameter :: SWITCH_FLOOR = 1.0e-8_dp
      !! pyscf zeroes a scaled distance below this before switching it

   type :: pcm_surface_t
      !! The discretised cavity: one Gaussian per surviving Lebedev point
      integer :: n_points = 0
      integer :: n_discarded = 0                !! buried points dropped
      integer, allocatable :: parent(:)         !! owning atom
      real(dp), allocatable :: coords(:, :)     !! (3, n), bohr
      real(dp), allocatable :: normals(:, :)    !! (3, n), outward unit normals
      real(dp), allocatable :: zeta(:)          !! the Gaussian's exponent is zeta**2
      real(dp), allocatable :: weights(:)       !! Lebedev weight times 4 pi
      real(dp), allocatable :: switch_f(:)      !! SWIG switching value, in (0, 1]
      real(dp), allocatable :: area(:)          !! weight * R**2 * switch, bohr**2
      real(dp), allocatable :: r_parent(:)      !! cavity radius of the owning atom
   end type pcm_surface_t

   type :: trc_pcm_t
      !! Everything an SCF iteration needs from the continuum
      integer :: model = PCM_MODEL_NONE
      real(dp) :: eps = 0.0_dp                  !! solvent dielectric constant
      real(dp) :: f_eps = 0.0_dp                !! the model's dielectric scale factor
      type(pcm_surface_t) :: surface
      real(dp), allocatable :: kmat(:, :)       !! K, LU-factored in place
      integer, allocatable :: ipiv(:)           !! its pivots
      real(dp), allocatable :: rmat(:, :)       !! R, IEF-PCM only
      real(dp), allocatable :: v_nuc(:)         !! nuclear potential at the points
      real(dp), allocatable :: v(:)             !! last solute potential, v_nuc - v_elec
      real(dp), allocatable :: q(:)             !! last symmetrised charges
      real(dp), allocatable :: wq(:)            !! -q: what the Fock kernel is handed
      real(dp) :: e_diel = 0.0_dp               !! q . v / 2 at the last density
      real(dp) :: q_total = 0.0_dp              !! total apparent surface charge
      type(trc_linalg_t), pointer :: la => null()
      !! the LU, solves and products: the caller's handles (the SCF's, made
      !! once per calculation) when init is given them, else PCM's own
      logical :: own_la = .false.               !! `la` was allocated here
      logical :: on_device = .false.            !! `this` is mapped; every array
                                                !! above follows it as it is made
   contains
      procedure :: init => pcm_init
      procedure :: build => pcm_build
      procedure :: build_operators => pcm_build_operators
      procedure :: solve => pcm_solve
      procedure :: operator => pcm_operator
      procedure :: release => pcm_release
   end type trc_pcm_t

contains

   pure function pcm_switch_h(x) result(h)
      !! The elementary switching function, eq. 3.19 of the SWIG paper
      !!
      !! A quintic rising from 0 to 1 across [0, 1] with two vanishing
      !! derivatives at each end. (The paper prints it with a typo; this is
      !! the corrected form pyscf uses.)
      real(dp), intent(in) :: x
      real(dp) :: h
      if (x <= 0.0_dp) then
         h = 0.0_dp
      else if (x >= 1.0_dp) then
         h = 1.0_dp
      else
         h = x**3*(10.0_dp - 15.0_dp*x + 6.0_dp*x**2)
      end if
   end function pcm_switch_h

   integer function pcm_model_from_name(name) result(model)
      !! PCM_MODEL_CPCM / PCM_MODEL_IEFPCM from a name, PCM_MODEL_NONE for
      !! an empty one, -1 for anything else
      character(len=*), intent(in) :: name
      select case (trim(adjustl(name)))
      case ("")
         model = PCM_MODEL_NONE
      case ("cpcm", "c-pcm", "CPCM", "C-PCM")
         model = PCM_MODEL_CPCM
      case ("iefpcm", "ief-pcm", "IEFPCM", "IEF-PCM")
         model = PCM_MODEL_IEFPCM
      case default
         model = -1
      end select
   end function pcm_model_from_name

   subroutine build_pcm_surface(atomic_numbers, coords, angular_points, radii_scale, surface, error)
      !! Tessellate the cavity: Lebedev points per atom, switched and pruned
      integer, intent(in) :: atomic_numbers(:)
      real(dp), intent(in) :: coords(:, :)        !! (3, natm), bohr
      integer, intent(in) :: angular_points       !! Lebedev points per sphere
      real(dp), intent(in) :: radii_scale         !! van der Waals radius to cavity radius
      type(pcm_surface_t), intent(out) :: surface
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: unit_points(:, :), unit_weights(:)
      real(dp), allocatable :: radii(:), r_switch(:), r_inner(:)
      real(dp), allocatable :: pt(:, :), nrm(:, :), zet(:), wgt(:), swf(:)
      integer, allocatable :: own(:)
      real(dp) :: point(3), xi, ratio, alpha, w, f, d
      integer :: natm, iatom, jatom, ipt, n_kept, ix, i
      character(len=16) :: num

      if (error%has_error()) return
      natm = size(atomic_numbers)
      if (size(coords, 1) /= 3 .or. size(coords, 2) /= natm) then
         call error%set(ERROR_VALIDATION, "PCM surface: coordinates must be (3, natm)")
         return
      end if

      ! The exponent table is fitted per Lebedev order; see PCM_XI_ORDERS.
      ix = 0
      do i = 1, N_XI
         if (PCM_XI_ORDERS(i) == angular_points) ix = i
      end do
      if (ix == 0 .or. .not. lebedev_is_available(angular_points)) then
         write (num, '(i0)') angular_points
         call error%set(ERROR_VALIDATION, "PCM: "//trim(num)//" points per sphere has no "// &
                        "fitted SWIG exponent [J. Chem. Phys. 122, 194110 (2005), Table II]; "// &
                        "use one of 6, 14, 26, 38, 50, 86, 110, 146, 170, 194, 302 (the "// &
                        "default), 350, 434, 590, 770, 974, ... 5810")
         return
      end if
      xi = XI_VALUES(ix)

      call lebedev_grid(angular_points, unit_points, unit_weights, error)
      if (error%has_error()) return

      ! The switching region of each sphere, eqs. 3.11-3.13 of the SWIG
      ! paper: a shell of width R_sw just inside the cavity radius, sized so
      ! that the region a point sweeps as it fades matches the point spacing.
      allocate (radii(natm), r_switch(natm), r_inner(natm))
      do iatom = 1, natm
         call pcm_cavity_radius(atomic_numbers(iatom), radii_scale, radii(iatom), error)
         if (error%has_error()) return
         r_switch(iatom) = radii(iatom)*sqrt(14.0_dp/real(angular_points, dp))
         ratio = radii(iatom)/r_switch(iatom)
         alpha = 0.5_dp + ratio - sqrt(ratio**2 - 1.0_dp/28.0_dp)
         r_inner(iatom) = radii(iatom) - alpha*r_switch(iatom)
      end do

      allocate (pt(3, natm*angular_points), nrm(3, natm*angular_points))
      allocate (zet(natm*angular_points), wgt(natm*angular_points))
      allocate (swf(natm*angular_points), own(natm*angular_points))

      n_kept = 0
      do iatom = 1, natm
         do ipt = 1, angular_points
            point = radii(iatom)*unit_points(:, ipt) + coords(:, iatom)
            w = 4.0_dp*PI*unit_weights(ipt)
            ! The product of every other sphere's elementary switch. One
            ! factor at zero -- the point is inside that sphere's inner
            ! radius -- kills the product, which is what "buried" means.
            f = 1.0_dp
            do jatom = 1, natm
               if (jatom == iatom) cycle
               d = (norm2(point - coords(:, jatom)) - r_inner(jatom))/r_switch(jatom)
               if (d < SWITCH_FLOOR) d = 0.0_dp
               f = f*pcm_switch_h(d)
            end do
            if (w*f <= BURIED_CUTOFF) then
               surface%n_discarded = surface%n_discarded + 1
               cycle
            end if
            n_kept = n_kept + 1
            pt(:, n_kept) = point
            nrm(:, n_kept) = unit_points(:, ipt)
            wgt(n_kept) = w
            swf(n_kept) = f
            own(n_kept) = iatom
            zet(n_kept) = xi/(radii(iatom)*sqrt(w))
         end do
      end do

      if (n_kept == 0) then
         call error%set(ERROR_VALIDATION, "PCM surface: every point was buried, which no "// &
                        "molecular cavity produces -- check that the geometry is in bohr")
         return
      end if

      surface%n_points = n_kept
      surface%parent = own(1:n_kept)
      surface%coords = pt(:, 1:n_kept)
      surface%normals = nrm(:, 1:n_kept)
      surface%zeta = zet(1:n_kept)
      surface%weights = wgt(1:n_kept)
      surface%switch_f = swf(1:n_kept)
      allocate (surface%area(n_kept), surface%r_parent(n_kept))
      do ipt = 1, n_kept
         surface%r_parent(ipt) = radii(own(ipt))
         surface%area(ipt) = wgt(ipt)*radii(own(ipt))**2*swf(ipt)
      end do
   end subroutine build_pcm_surface

   pure real(dp) function s_elem(i, j, n, coords, zeta, switch_f)
      !! One element of S: erf(zeta_ij r)/r off the diagonal, zeta sqrt(2/pi)/F on it
      !$acc routine seq
      integer, intent(in) :: i, j, n
      real(dp), intent(in) :: coords(3, n), zeta(n), switch_f(n)
      real(dp) :: dx, dy, dz, r, zij
      if (i == j) then
         s_elem = zeta(i)*sqrt(2.0_dp/PI)/switch_f(i)
      else
         dx = coords(1, i) - coords(1, j)
         dy = coords(2, i) - coords(2, j)
         dz = coords(3, i) - coords(3, j)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         zij = zeta(i)*zeta(j)/sqrt(zeta(i)**2 + zeta(j)**2)
         s_elem = erf(zij*r)/r
      end if
   end function s_elem

   pure real(dp) function d_elem(i, j, n, coords, normals, zeta, r_parent)
      !! One element of D, S's derivative along the normal of point j
      !$acc routine seq
      integer, intent(in) :: i, j, n
      real(dp), intent(in) :: coords(3, n), normals(3, n), zeta(n), r_parent(n)
      real(dp) :: dx, dy, dz, r, zij, xr, nr
      if (i == j) then
         d_elem = -zeta(i)*sqrt(2.0_dp/PI)/(2.0_dp*r_parent(i))
      else
         dx = coords(1, i) - coords(1, j)
         dy = coords(2, i) - coords(2, j)
         dz = coords(3, i) - coords(3, j)
         r = sqrt(dx*dx + dy*dy + dz*dz)
         zij = zeta(i)*zeta(j)/sqrt(zeta(i)**2 + zeta(j)**2)
         xr = zij*r
         nr = dx*normals(1, j) + dy*normals(2, j) + dz*normals(3, j)
         d_elem = erf(xr)/r*nr/r**2 - 2.0_dp*xr/sqrt(PI)*exp(-xr*xr)*nr/r**3
      end if
   end function d_elem

   subroutine s_fill(n, coords, zeta, switch_f, smat)
      !! S, on arrays already where the build computes (the device, on GPU)
      !!
      !! Every element is independent, so the matrix is one `do concurrent`
      !! over the square: at 10^4 points that is 10^8 erf evaluations.
      integer, intent(in) :: n
      real(dp), intent(in) :: coords(3, n), zeta(n), switch_f(n)
      real(dp), intent(out) :: smat(n, n)
      integer :: i, j
      do concurrent(j = 1:n, i = 1:n)
         smat(i, j) = s_elem(i, j, n, coords, zeta, switch_f)
      end do
   end subroutine s_fill

   subroutine da_fill(n, coords, normals, zeta, r_parent, area, dmat)
      !! D A: D with its columns scaled by the point areas, as s_fill
      integer, intent(in) :: n
      real(dp), intent(in) :: coords(3, n), normals(3, n), zeta(n), r_parent(n), area(n)
      real(dp), intent(out) :: dmat(n, n)
      integer :: i, j
      do concurrent(j = 1:n, i = 1:n)
         dmat(i, j) = d_elem(i, j, n, coords, normals, zeta, r_parent)*area(j)
      end do
   end subroutine da_fill

   subroutine pcm_sd_matrices(n, coords, normals, zeta, switch_f, r_parent, smat, dmat)
      !! S, and D when it is present, for a surface; host arrays in and out
      !!
      !! The element builders K is formed from (s_fill, da_fill), wrapped for
      !! a caller that wants the matrices themselves -- check_pcm_cavity
      !! compares them with pyscf's. D here is D itself: unit areas.
      integer, intent(in) :: n
      real(dp), intent(in) :: coords(3, n), normals(3, n), zeta(n), switch_f(n), r_parent(n)
      real(dp), intent(out) :: smat(n, n)
      real(dp), intent(out), optional :: dmat(n, n)
      real(dp), allocatable :: ones(:)
      !$acc enter data copyin(coords, normals, zeta, switch_f, r_parent) create(smat)
      call s_fill(n, coords, zeta, switch_f, smat)
      !$acc update self(smat)
      !$acc exit data delete(smat)
      if (present(dmat)) then
         allocate (ones(n))
         ones = 1.0_dp
         !$acc enter data copyin(ones) create(dmat)
         call da_fill(n, coords, normals, zeta, r_parent, ones, dmat)
         !$acc update self(dmat)
         !$acc exit data delete(ones, dmat)
      end if
      !$acc exit data delete(coords, normals, zeta, switch_f, r_parent)
   end subroutine pcm_sd_matrices

   subroutine pcm_init(this, model, eps, error, la)
      !! Settle the model and its dielectric scale factor
      class(trc_pcm_t), intent(inout) :: this
      integer, intent(in) :: model           !! PCM_MODEL_CPCM or PCM_MODEL_IEFPCM
      real(dp), intent(in) :: eps            !! solvent dielectric constant
      type(error_t), intent(inout) :: error
      type(trc_linalg_t), intent(inout), target, optional :: la
      !! cuBLAS/cuSOLVER handles to borrow (trc_linalg's init is create-once,
      !! so PCM adds none of its own); absent, PCM makes and releases its own

      if (error%has_error()) return
      call this%release()
      ! Required rather than defaulted: a default would silently pick a solvent.
      if (eps <= 1.0_dp) then
         call error%set(ERROR_VALIDATION, "PCM: the solvent dielectric constant must be "// &
                        "greater than one (water is 78.3553)")
         return
      end if
      select case (model)
      case (PCM_MODEL_CPCM)
         this%f_eps = (eps - 1.0_dp)/eps
      case (PCM_MODEL_IEFPCM)
         this%f_eps = (eps - 1.0_dp)/(eps + 1.0_dp)
      case default
         call error%set(ERROR_VALIDATION, "PCM: the model must be 'cpcm' (conductor-like) "// &
                        "or 'iefpcm' (integral equation formalism); they are different "// &
                        "models with different energies, so neither stands in for the other")
         return
      end select
      this%eps = eps
      this%model = model
      ! Parent first, as the basis does (trc_api's basis_to_device); each
      ! array is then mapped as it is made, and released in reverse. The
      ! handles live outside `this` (a pointer, to the caller's object or to
      ! one allocated here), so their device scalars are never part of
      ! `this`'s mapping.
      !$acc enter data copyin(this)
      if (present(la)) then
         this%la => la
         this%own_la = .false.
      else
         allocate (this%la)
         this%own_la = .true.
      end if
      call this%la%init(1)
      this%on_device = .true.
   end subroutine pcm_init

   subroutine pcm_build(this, b, angular_points, radii_scale, error)
      !! Cavity, operators and nuclear potential, once per geometry
      class(trc_pcm_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: angular_points
      real(dp), intent(in) :: radii_scale
      type(error_t), intent(inout) :: error

      integer :: k, ia, n
      real(dp) :: d

      if (error%has_error()) return
      if (this%model == PCM_MODEL_NONE) then
         call error%set(ERROR_VALIDATION, "PCM: build before init"); return
      end if
      if (maxval(b%sh_l) > PCM_LMAX) then
         call error%set(ERROR_VALIDATION, "PCM: the basis has angular momentum above this "// &
                        "build's PCM surface-integral kernels (gen_1e.py --pcm --lmax)")
         return
      end if
      call build_pcm_surface(nint(b%at_z), b%at_r, angular_points, radii_scale, this%surface, error)
      if (error%has_error()) return
      call this%build_operators(error)
      if (error%has_error()) return

      ! The nuclear potential each smeared point feels. A normalised
      ! Gaussian of exponent zeta**2 against a point nucleus integrates to
      ! erf(zeta d)/d exactly.
      n = this%surface%n_points
      allocate (this%v_nuc(n))
      this%v_nuc = 0.0_dp
      do k = 1, n
         do ia = 1, b%natm
            d = norm2(this%surface%coords(:, k) - b%at_r(:, ia))
            this%v_nuc(k) = this%v_nuc(k) + b%at_z(ia)*erf(this%surface%zeta(k)*d)/d
         end do
      end do
      !$acc enter data copyin(this%v_nuc)
   end subroutine pcm_build

   subroutine pcm_build_operators(this, error)
      !! Form K from the surface, factor it, and keep R where the model has one
      !!
      !! On the device under the GPU build, from the S and D builds through
      !! the D A S product to the LU; only the pivot status comes back.
      class(trc_pcm_t), intent(inout) :: this
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: smat(:, :)
      integer :: n, info
      real(dp) :: f, f2pi
      character(len=16) :: num

      if (error%has_error()) return
      if (.not. this%on_device) then
         call error%set(ERROR_VALIDATION, "PCM: build_operators before init"); return
      end if
      n = this%surface%n_points
      f = this%f_eps
      f2pi = f/(2.0_dp*PI)
      allocate (this%kmat(n, n), this%ipiv(n), this%v(n), this%q(n), this%wq(n))
      this%v = 0.0_dp; this%q = 0.0_dp; this%wq = 0.0_dp
      !$acc enter data copyin(this%surface%coords, this%surface%normals, this%surface%zeta) &
      !$acc            copyin(this%surface%switch_f, this%surface%r_parent, this%surface%area) &
      !$acc            create(this%kmat, this%ipiv) copyin(this%v, this%q, this%wq)
      associate (sf => this%surface)
         if (this%model == PCM_MODEL_IEFPCM) then
            ! S into a scratch copy, D A into R's storage; then
            ! K = S - f/(2 pi) (D A) S, and R = -f (I - D A/(2 pi)) in place.
            allocate (smat(n, n), this%rmat(n, n))
            !$acc enter data create(smat, this%rmat)
            call s_fill(n, sf%coords, sf%zeta, sf%switch_f, smat)
            call da_fill(n, sf%coords, sf%normals, sf%zeta, sf%r_parent, sf%area, this%rmat)
            call copy_square(n, smat, this%kmat)
            call this%la%gemm('N', 'N', n, n, n, -f2pi, this%rmat, n, smat, n, 1.0_dp, this%kmat, n)
            !$acc exit data delete(smat)
            deallocate (smat)
            call scale_shift(n, this%rmat, f2pi, -f)
         else
            call s_fill(n, sf%coords, sf%zeta, sf%switch_f, this%kmat)
         end if
      end associate
      call this%la%getrf(n, this%kmat, n, this%ipiv, info)
      if (info /= 0) then
         write (num, '(i0)') info
         call error%set(ERROR_VALIDATION, "PCM: K is singular (getrf info = "//trim(num)// &
                        "), which a physical cavity does not produce; two surface points "// &
                        "may coincide -- check the geometry")
      end if
   end subroutine pcm_build_operators

   subroutine copy_square(n, a, b)
      integer, intent(in) :: n
      real(dp), intent(in) :: a(n, n)
      real(dp), intent(out) :: b(n, n)
      integer :: i, j
      do concurrent(j = 1:n, i = 1:n)
         b(i, j) = a(i, j)
      end do
   end subroutine copy_square

   subroutine axpby(n, alpha, x, beta, y)
      !! y := alpha x + beta y, where the arrays are (the device, on GPU)
      !!
      !! beta = 0 does not read y, as BLAS does not: y may be freshly
      !! allocated, and 0 times whatever NaN the heap left there is NaN.
      !! (ifx found that one: IEF-PCM's scratch vector, a dozen SCFs into
      !! a process, and the eigensolver refusing the Fock matrix after.)
      integer, intent(in) :: n
      real(dp), intent(in) :: alpha, beta
      real(dp), intent(in) :: x(n)
      real(dp), intent(inout) :: y(n)
      integer :: k
      if (beta == 0.0_dp) then
         do concurrent(k = 1:n)
            y(k) = alpha*x(k)
         end do
      else
         do concurrent(k = 1:n)
            y(k) = alpha*x(k) + beta*y(k)
         end do
      end if
   end subroutine axpby

   subroutine scale_shift(n, a, alpha, beta)
      !! a := alpha a + beta I
      integer, intent(in) :: n
      real(dp), intent(inout) :: a(n, n)
      real(dp), intent(in) :: alpha, beta
      integer :: i, j
      do concurrent(j = 1:n, i = 1:n)
         a(i, j) = alpha*a(i, j) + merge(beta, 0.0_dp, i == j)
      end do
   end subroutine scale_shift

   subroutine pcm_solve(this, v, q)
      !! Symmetrised apparent surface charges for one solute potential
      !!
      !! `K q = R v`, and for IEF-PCM's nonsymmetric K also `q_t = R^T K^-T v`,
      !! averaged as pyscf does. C-PCM's K is symmetric and R a scalar, so
      !! one solve is both. `v` and `q` live where K does (present on the
      !! device under the GPU build). The one place K is used per iteration:
      !! a matrix-free solve replaces this and build_operators, nothing else.
      class(trc_pcm_t), intent(inout) :: this
      real(dp), intent(in) :: v(this%surface%n_points)
      real(dp), intent(inout) :: q(this%surface%n_points)

      real(dp), allocatable :: adj(:)
      integer :: n
      real(dp) :: f

      n = this%surface%n_points
      f = this%f_eps
      if (this%model == PCM_MODEL_IEFPCM) then
         call this%la%gemv('N', n, n, 1.0_dp, this%rmat, n, v, 0.0_dp, q)
         call this%la%getrs('N', n, 1, this%kmat, n, this%ipiv, q, n)
         allocate (adj(n))
         !$acc enter data create(adj)
         call axpby(n, 1.0_dp, v, 0.0_dp, adj)
         call this%la%getrs('T', n, 1, this%kmat, n, this%ipiv, adj, n)
         ! q := (q + R^T adj)/2
         call this%la%gemv('T', n, n, 0.5_dp, this%rmat, n, adj, 0.5_dp, q)
         !$acc exit data delete(adj)
      else
         call axpby(n, -f, v, 0.0_dp, q)
         call this%la%getrs('N', n, 1, this%kmat, n, this%ipiv, q, n)
      end if
   end subroutine pcm_solve

   subroutine pcm_point_range(n, rank, nranks, k0, k1)
      !! Rank `rank` (0-based) of `nranks`'s contiguous share of points 1..n
      !!
      !! Balanced to within one point: the first `mod(n, nranks)` ranks get
      !! `n/nranks + 1`, the rest `n/nranks`. More ranks than points leaves
      !! some ranks with an empty share, signalled by k1 < k0 -- the
      !! caller's job to skip, not to call the integral drivers with zero
      !! points.
      integer, intent(in) :: n, rank, nranks
      integer, intent(out) :: k0, k1
      integer :: base, extra
      base = n/nranks
      extra = mod(n, nranks)
      if (rank < extra) then
         k0 = rank*(base + 1) + 1
         k1 = k0 + base
      else
         k0 = extra*(base + 1) + (rank - extra)*base + 1
         k1 = k0 + base - 1
      end if
   end subroutine pcm_point_range

   subroutine pcm_sum_ranks(comm, g, n)
      !! Sum a contiguous array over the ranks, in place, then rank 0's
      !! copy everywhere -- an explicit-shape `g(n)` so a 2-D actual (vmat)
      !! and a 1-D one (this%v) both bind by element sequence, the same
      !! reduce-then-broadcast trc_scf's sum_flat and trc_fock's sum_ranks
      !! use for the other matrices an iteration reduces.
      type(comm_t), intent(in) :: comm
      integer, intent(in) :: n
      real(dp), intent(inout) :: g(n)
      call allreduce(comm, g, op=MPI_SUM)
      call bcast(comm, g, n, 0)
   end subroutine pcm_sum_ranks

   subroutine pcm_operator(this, b, pl, dtot, vmat, e_diel, comm)
      !! One SCF iteration's continuum: charges from this density, then their operator
      !!
      !! `dtot` (the total density, alpha plus beta) and `vmat` are Cartesian
      !! and device-resident. `vmat` is overwritten with the full operator
      !! the charges add to the Fock matrix; `e_diel` already carries the
      !! half. Everything stays on the device; the charges come back (n
      !! numbers) only for their total.
      !!
      !! `comm`, when present with more than one rank, splits the two
      !! integral passes -- trc_pcm_potential and trc_pcm_fock, the only
      !! parts of an iteration that scale with the integral work rather
      !! than with the (replicated, cheap) cavity solve -- over a
      !! contiguous range of surface points per rank (pcm_point_range).
      !! Each rank computes its own range into the full-length vector,
      !! zero elsewhere; summing those across ranks (adding exact zeros)
      !! reproduces the single-rank vector, and the broadcast after makes
      !! every rank's bits agree, the same reduce-then-broadcast
      !! trc_fock's sum_ranks uses for the real Fock matrix. The solve, the
      !! energy and everything after need no communication: once `this%v`
      !! agrees, they are identical, replicated work on every rank.
      !! Absent `comm`, or a single rank, this is exactly the unsplit path.
      class(trc_pcm_t), intent(inout) :: this
      type(trc_basis_t), intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      real(dp), intent(in) :: dtot(b%nao, b%nao)
      real(dp), intent(out) :: vmat(b%nao, b%nao)
      real(dp), intent(out) :: e_diel
      type(comm_t), intent(in), optional :: comm
      integer :: n, i, j, k0, k1, nranks

      n = this%surface%n_points
      nranks = 1
      if (present(comm)) nranks = comm%size()

      if (nranks > 1) then
         call pcm_point_range(n, comm%rank(), nranks, k0, k1)
         do concurrent(i = 1:n)
            this%v(i) = 0.0_dp
         end do
         if (k1 >= k0) call trc_pcm_potential(b, pl, k1 - k0 + 1, this%surface%coords(:, k0:k1), &
                                              this%surface%zeta(k0:k1), dtot, this%v(k0:k1))
         !$acc update self(this%v)
         call pcm_sum_ranks(comm, this%v, n)
         !$acc update device(this%v)
      else
         call trc_pcm_potential(b, pl, n, this%surface%coords, this%surface%zeta, dtot, this%v)
      end if
      ! Nuclei push the potential up, electrons pull it down.
      call axpby(n, 1.0_dp, this%v_nuc, -1.0_dp, this%v)
      call this%solve(this%v, this%q)
      this%e_diel = 0.5_dp*this%la%dot(n, this%q, this%v)
      e_diel = this%e_diel
      ! The electron's charge against the apparent charges' potential.
      call axpby(n, -1.0_dp, this%q, 0.0_dp, this%wq)
      do concurrent(i = 1:b%nao, j = 1:b%nao)
         vmat(i, j) = 0.0_dp
      end do
      if (nranks > 1) then
         ! Same range as the potential pass, above; each rank's call adds
         ! only its share of points into vmat (zeroed just above it).
         if (k1 >= k0) call trc_pcm_fock(b, pl, k1 - k0 + 1, this%surface%coords(:, k0:k1), &
                                         this%surface%zeta(k0:k1), this%wq(k0:k1), vmat)
         !$acc update self(vmat)
         call pcm_sum_ranks(comm, vmat, b%nao*b%nao)
         !$acc update device(vmat)
      else
         call trc_pcm_fock(b, pl, n, this%surface%coords, this%surface%zeta, this%wq, vmat)
      end if
      !$acc update self(this%q)
      this%q_total = sum(this%q)
   end subroutine pcm_operator

   subroutine pcm_release(this)
      !! Release everything; reusable after another `init`
      class(trc_pcm_t), intent(inout) :: this
      if (this%on_device) then
         if (allocated(this%kmat)) then
            !$acc exit data delete(this%surface%coords, this%surface%normals, this%surface%zeta) &
            !$acc           delete(this%surface%switch_f, this%surface%r_parent, this%surface%area) &
            !$acc           delete(this%kmat, this%ipiv, this%v, this%q, this%wq)
         end if
         if (allocated(this%rmat)) then
            !$acc exit data delete(this%rmat)
         end if
         if (allocated(this%v_nuc)) then
            !$acc exit data delete(this%v_nuc)
         end if
         if (this%own_la) then
            call this%la%release()
            deallocate (this%la)
         end if
         nullify (this%la)
         this%own_la = .false.
         !$acc exit data delete(this)
         this%on_device = .false.
      end if
      if (allocated(this%kmat)) deallocate (this%kmat)
      if (allocated(this%ipiv)) deallocate (this%ipiv)
      if (allocated(this%rmat)) deallocate (this%rmat)
      if (allocated(this%v_nuc)) deallocate (this%v_nuc)
      if (allocated(this%v)) deallocate (this%v)
      if (allocated(this%q)) deallocate (this%q)
      if (allocated(this%wq)) deallocate (this%wq)
      this%surface = pcm_surface_t()
      this%model = PCM_MODEL_NONE
      this%eps = 0.0_dp; this%f_eps = 0.0_dp
      this%e_diel = 0.0_dp; this%q_total = 0.0_dp
   end subroutine pcm_release

end module trc_pcm
