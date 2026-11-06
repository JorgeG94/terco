!
! CHELPG atomic charges: point charges fitted to the molecule's own
! electrostatic potential on a grid outside it
!
! Ported from metalquicha's mqc_czt_charges (chelpg_grid, chelpg_charges),
! mqc_ -> trc_ renamed. The grid rule, the radii and the constrained fit are
! mqc's; the potential is terco's own PCM surface-integral kernel with the
! Gaussian smearing taken to the point-charge limit, and the grid and the
! fit's normal equations are `do concurrent` kernels.
!
module trc_chelpg
   !! Charges `q_A` on the nuclei that best reproduce the molecule's
   !! electrostatic potential on a shell of points outside it:
   !!
   !!     minimise  sum_g ( V(g) - sum_A q_A / |r_g - R_A| )^2
   !!     subject to  sum_A q_A = Q
   !!
   !! with `V(g) = sum_A Z_A/|r_g - R_A| - sum_{mu nu} D_mu,nu <mu|1/|r - r_g||nu>`
   !! (Breneman and Wiberg, J. Comput. Chem. 11, 361 (1990)). The potential
   !! is an observable of the total density, so the charges are much less
   !! basis-set dependent than a population analysis, and they are the ones
   !! to embed a fragment with: reproducing a field is what they are fitted
   !! to do. R and U references alike: only the total density enters.
   !!
   !! **The grid** is Breneman and Wiberg's, exactly as mqc builds it: a cubic
   !! lattice of spacing `CHELPG_SPACING` (0.3 Angstrom) over the molecule's
   !! bounding box grown by `CHELPG_HEAD_SPACE` (2.8 Angstrom) on every side,
   !! starting at the box's low corner, keeping a point only when it lies
   !! outside every atom's van der Waals sphere and within the head space of
   !! at least one NUCLEUS. Inside a sphere a point charge's potential
   !! diverges while the real one does not; far outside, every charge
   !! distribution with the right multipoles looks alike.
   !!
   !! This is NOT GPU4PySCF's `eval_chelpg_layer_gpu` grid: that one (after
   !! Q-Chem) measures the head space from the van der Waals SURFACE rather
   !! than from the nucleus, grows the box by the head space plus the
   !! outermost atom's radius, and by default smooths both cutoffs
   !! (`ifqchem`). Same method, different points, different charges: on
   !! water/6-31G at mqc's geometry, 2.7e-4 e from GPU4PySCF's default and
   !! 5.7e-4 e from its unsmoothed variant (measured, GPU4PySCF 1.9.0, pyscf
   !! 2.14). terco follows mqc.
   !!
   !! **The fit** carries the total-charge constraint only, through one
   !! Lagrange multiplier, as mqc's does: no dipole constraint. The
   !! `natm+1` system is symmetric but indefinite, so it is solved by LU.
   !!
   !! Refused: a ghost atom (charge zero -- whether it should carry a fitted
   !! charge, and whether its sphere should exclude points, is a convention
   !! terco has not chosen), an element without a tabulated radius, and a
   !! basis above the potential kernel's `PCM_LMAX`. terco reads
   !! all-electron bases only (the JSON reader refuses an ECP), so there is no
   !! core potential to account for.
   use pic_types, only: dp, default_int
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use pic_lapack_interfaces, only: pic_getrf, pic_getrs
   use trc_constants, only: ANGSTROM_TO_BOHR
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_pcm_potential, PCM_LMAX
   implicit none
   private

   public :: trc_chelpg_charges
   public :: trc_chelpg_grid
   public :: trc_chelpg_potential
   public :: trc_chelpg_vdw_radius
   public :: trc_chelpg_options_t, trc_chelpg_fit_t
   public :: CHELPG_SPACING, CHELPG_HEAD_SPACE

   real(dp), parameter :: CHELPG_SPACING = 0.3_dp
      !! Lattice spacing, Angstrom: Breneman and Wiberg's value, and mqc's.
   real(dp), parameter :: CHELPG_HEAD_SPACE = 2.8_dp
      !! Angstrom. Both the box margin and the cutoff: the lattice extends
      !! this far beyond the outermost nuclei, and a point further than this
      !! from every nucleus is dropped. Breneman and Wiberg's value, and mqc's.

   real(dp), parameter :: POINT_ZETA = 1.0e150_dp
      !! The PCM kernel's operator is erf(zeta r)/r. At this zeta its
      !! attenuation theta = zeta^2/(zeta^2 + p) is 1.0 to the last bit for
      !! every primitive pair exponent p below ~1e284, so the kernel returns
      !! the plain 1/r nuclear-attraction integral (its own "theta -> 1 is
      !! the nuclear attraction exactly"), and zeta^2 = 1e300 stays finite.

   integer, parameter :: N_VDW = 96
   real(dp), parameter :: VDW_RADII_ANGSTROM(N_VDW) = [ &
                          1.10_dp, 1.40_dp, 1.81_dp, 1.53_dp, 1.92_dp, 1.70_dp, 1.55_dp, 1.52_dp, &
                          1.47_dp, 1.54_dp, 2.27_dp, 1.73_dp, 1.84_dp, 2.10_dp, 1.80_dp, 1.80_dp, &
                          1.75_dp, 1.88_dp, 2.75_dp, 2.31_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, &
                          0.00_dp, 0.00_dp, 0.00_dp, 1.63_dp, 1.40_dp, 1.39_dp, 1.87_dp, 2.11_dp, &
                          1.85_dp, 1.90_dp, 1.85_dp, 2.02_dp, 3.03_dp, 2.49_dp, 0.00_dp, 0.00_dp, &
                          0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 1.63_dp, 1.72_dp, 1.58_dp, &
                          1.93_dp, 2.17_dp, 2.06_dp, 2.06_dp, 1.98_dp, 2.16_dp, 3.43_dp, 2.68_dp, &
                          0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, &
                          0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, &
                          0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp, 1.75_dp, 1.66_dp, 1.55_dp, &
                          1.96_dp, 2.02_dp, 2.07_dp, 1.97_dp, 2.02_dp, 2.20_dp, 3.48_dp, 2.83_dp, &
                          0.00_dp, 0.00_dp, 0.00_dp, 1.86_dp, 0.00_dp, 0.00_dp, 0.00_dp, 0.00_dp]
      !! Van der Waals radii, Angstrom, indexed by atomic number: mqc's
      !! `mqc_elements` table verbatim. Bondi [J. Phys. Chem. 68, 441 (1964)],
      !! with Rowland and Taylor's 1.10 for hydrogen [J. Phys. Chem. 100,
      !! 7384 (1996)] and Mantina et al. [J. Phys. Chem. A 113, 5806 (2009)]
      !! for the main-group elements Bondi did not cover (and Li, 1.81). Zero
      !! where neither tabulates one; such an element is refused. These are a
      !! convention the charges depend on: a different set moves the points
      !! the fit sees. (PCM's own table, trc_pcm_radii, is pyscf's modified
      !! Bondi: the same H-Ca except Li at 1.82.)

   type :: trc_chelpg_options_t
      !! The grid's two lengths, Angstrom; the defaults are Breneman and Wiberg's.
      real(dp) :: spacing = CHELPG_SPACING
      real(dp) :: head_space = CHELPG_HEAD_SPACE
   end type trc_chelpg_options_t

   type :: trc_chelpg_fit_t
      !! How well the fitted charges reproduce the potential they were fitted to
      integer :: n_points = 0      !! grid points the fit used
      real(dp) :: rms = 0.0_dp     !! rms of V - V_charges over the grid, Eh/e
      real(dp) :: rrms = 0.0_dp
         !! `rms` relative to the rms of V itself (the residual a zero-charge
         !! model would leave): the usual measure of an ESP fit, below ~0.1
         !! for a polar molecule.
   end type trc_chelpg_fit_t

contains

   pure function trc_chelpg_vdw_radius(atomic_number) result(radius)
      !! The CHELPG exclusion radius of an element, Bohr; 0 where none is tabulated
      integer, intent(in) :: atomic_number
      real(dp) :: radius
      radius = 0.0_dp
      if (atomic_number >= 1 .and. atomic_number <= N_VDW) then
         radius = VDW_RADII_ANGSTROM(atomic_number)*ANGSTROM_TO_BOHR
      end if
   end function trc_chelpg_vdw_radius

   subroutine trc_chelpg_charges(b, density, total_charge, charges, error, options, fit)
      !! The CHELPG charges of the density `density`, one per atom of `b`
      !!
      !! `density` is the total, both spins, in terco's CARTESIAN AO space
      !! (b%nao), the space `trc_hirshfeld_populations` reads too: a caller
      !! whose SCF ran spherical converts with `trc_sph`'s `sph_sph_to_cart`
      !! first. `total_charge` is the molecular charge the charges are
      !! constrained to sum to. `b` may be host-only or on the device.
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: density(b%nao, b%nao)
      real(dp), intent(in) :: total_charge
      real(dp), allocatable, intent(out) :: charges(:)   !! (natm), e
      type(error_t), intent(inout) :: error
      type(trc_chelpg_options_t), intent(in), optional :: options
      type(trc_chelpg_fit_t), intent(out), optional :: fit

      real(dp), allocatable :: points(:, :), potential(:), residual(:)
      real(dp), allocatable :: gram(:), projection(:), normal(:, :), rhs(:, :)
      integer(default_int), allocatable :: pivots(:)
      integer(default_int) :: info
      real(dp) :: v_rms
      integer :: n, npt

      if (error%has_error()) return
      n = b%natm

      call trc_chelpg_grid(b, points, error, options)
      if (error%has_error()) return
      npt = size(points, 2)

      call trc_chelpg_potential(b, density, points, potential, error)
      if (error%has_error()) return

      ! The bordered normal equations: rows 1..n are sum_g (1/r_gA)(1/r_gB)
      ! and sum_g V(g)/r_gA, row and column n+1 the charge constraint.
      allocate (gram(n*n), projection(n))
      call chelpg_normal_kernel(n, b%at_r, npt, points, potential, gram, projection)
      allocate (normal(n + 1, n + 1), rhs(n + 1, 1), pivots(n + 1))
      normal(1:n, 1:n) = reshape(gram, [n, n])
      normal(1:n, n + 1) = 1.0_dp
      normal(n + 1, 1:n) = 1.0_dp
      normal(n + 1, n + 1) = 0.0_dp
      rhs(1:n, 1) = projection
      rhs(n + 1, 1) = total_charge

      call pic_getrf(normal, pivots, info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "chelpg charges: the fit matrix is singular, which "// &
                        "means two atoms are indistinguishable to every grid point")
         return
      end if
      call pic_getrs(normal, pivots, rhs, info=info)
      if (info /= 0) then
         call error%set(ERROR_VALIDATION, "chelpg charges: the constrained fit would not solve")
         return
      end if
      allocate (charges(n), source=rhs(1:n, 1))

      if (present(fit)) then
         allocate (residual(npt))
         call chelpg_residual_kernel(n, b%at_r, charges, npt, points, potential, residual)
         v_rms = sqrt(sum(potential**2)/real(npt, dp))
         fit%n_points = npt
         fit%rms = sqrt(sum(residual**2)/real(npt, dp))
         fit%rrms = 0.0_dp
         if (v_rms > 0.0_dp) fit%rrms = fit%rms/v_rms
      end if
      call logger%verbose("  chelpg: fitted "//to_string(n)//" charges to "// &
                          to_string(npt)//" grid points")
   end subroutine trc_chelpg_charges

   subroutine trc_chelpg_grid(b, points, error, options)
      !! The points a CHELPG fit is made on, Bohr, in lattice order (x fastest)
      !!
      !! mqc's `chelpg_grid` exactly: the lattice starts at the bounding
      !! box's low corner `min(R) - head_space` and runs
      !! `int((max - min + 2 head_space)/spacing) + 1` points per axis; a
      !! point is kept when no nucleus is closer than its van der Waals
      !! radius and the nearest one is no further than `head_space`.
      type(trc_basis_t), intent(in) :: b
      real(dp), allocatable, intent(out) :: points(:, :)   !! (3, n_points)
      type(error_t), intent(inout) :: error
      type(trc_chelpg_options_t), intent(in), optional :: options

      type(trc_chelpg_options_t) :: opts
      real(dp), allocatable :: radii(:)
      logical, allocatable :: kept(:)
      real(dp) :: low(3), high(3)
      real(dp) :: h, margin
      integer :: nx(3)
      integer :: ia, d, n_lattice, n_kept, g, i, j, k

      if (error%has_error()) return
      if (present(options)) opts = options
      if (opts%spacing <= 0.0_dp .or. opts%head_space <= 0.0_dp) then
         call error%set(ERROR_VALIDATION, "chelpg grid: the spacing and the head space must be positive")
         return
      end if
      h = opts%spacing*ANGSTROM_TO_BOHR
      margin = opts%head_space*ANGSTROM_TO_BOHR

      allocate (radii(b%natm))
      do ia = 1, b%natm
         if (b%at_z(ia) == 0.0_dp) then
            call error%set(ERROR_VALIDATION, "chelpg charges: atom "//to_string(ia)// &
                           " is a ghost (no nuclear charge); whether a ghost centre carries a "// &
                           "fitted charge and excludes grid points is not decided, so ghosts "// &
                           "are refused")
            return
         end if
         radii(ia) = trc_chelpg_vdw_radius(b%at_elem(ia))
         if (radii(ia) <= 0.0_dp) then
            call error%set(ERROR_VALIDATION, "chelpg charges: no van der Waals radius for element "// &
                           to_string(b%at_elem(ia))//", so the region excluded around it is undefined")
            return
         end if
      end do

      do d = 1, 3
         low(d) = minval(b%at_r(d, :)) - margin
         high(d) = maxval(b%at_r(d, :)) + margin
         nx(d) = int((high(d) - low(d))/h) + 1
      end do
      n_lattice = nx(1)*nx(2)*nx(3)

      allocate (kept(n_lattice))
      call chelpg_lattice_kernel(n_lattice, nx(1), nx(2), low, h, margin, b%natm, b%at_r, radii, kept)

      ! Compacted on the host, in lattice order, so the points (and the fit's
      ! sums over them) come out in mqc's own loop order.
      n_kept = count(kept)
      if (n_kept < b%natm + 1) then
         call error%set(ERROR_VALIDATION, "chelpg grid: kept "//to_string(n_kept)// &
                        " points, fewer than the "//to_string(b%natm + 1)//" unknowns the "// &
                        "fit has; use a finer spacing or a larger head space")
         return
      end if
      allocate (points(3, n_kept))
      n_kept = 0
      do g = 1, n_lattice
         if (.not. kept(g)) cycle
         call lattice_indices(g, nx(1), nx(2), i, j, k)
         n_kept = n_kept + 1
         points(1, n_kept) = low(1) + i*h
         points(2, n_kept) = low(2) + j*h
         points(3, n_kept) = low(3) + k*h
      end do
   end subroutine trc_chelpg_grid

   subroutine trc_chelpg_potential(b, density, points, potential, error)
      !! The molecule's electrostatic potential at `points`, Eh/e: nuclei minus electrons
      !!
      !! The electrons' part is `trc_pcm_potential`, PCM's own surface-
      !! integral kernel, at `POINT_ZETA` -- the point-charge limit of its
      !! erf(zeta r)/r operator. `density` as for `trc_chelpg_charges`.
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: density(b%nao, b%nao)
      real(dp), intent(in) :: points(:, :)            !! (3, n), Bohr
      real(dp), allocatable, intent(out) :: potential(:)
      type(error_t), intent(inout) :: error

      type(trc_basis_t) :: b_device
      real(dp), allocatable :: electronic(:)
      integer :: npt

      if (error%has_error()) return
      if (maxval(b%sh_l) > PCM_LMAX) then
         call error%set(ERROR_VALIDATION, "chelpg charges: the basis has l = "// &
                        to_string(maxval(b%sh_l))//", above this build's potential-integral "// &
                        "ceiling PCM_LMAX = "//to_string(PCM_LMAX)//" (TERCO_LMAX)")
         return
      end if
      npt = size(points, 2)
      allocate (electronic(npt), potential(npt))
      if (b%on_device) then
         call electronic_potential(b, npt, points, density, electronic)
      else
         ! A caller outside an SCF context may hand over a host-only basis;
         ! the kernel reads it on the device.
         b_device = b
         b_device%on_device = .false.
         call b_device%to_device()
         call electronic_potential(b_device, npt, points, density, electronic)
         call b_device%release()
      end if
      call chelpg_nuclear_kernel(b%natm, b%at_r, b%at_z, npt, points, electronic, potential)
   end subroutine trc_chelpg_potential

   subroutine electronic_potential(b, npt, pts, density, electronic)
      !! sum_{mu nu} D_mu,nu <mu|1/|r - r_g||nu> at every point; `b` on the device
      type(trc_basis_t), intent(in) :: b
      integer, intent(in) :: npt
      real(dp), intent(in) :: pts(3, npt), density(b%nao, b%nao)
      real(dp), intent(out) :: electronic(npt)

      type(trc_pairlist_t) :: pairs
      real(dp), allocatable :: zeta(:)

      allocate (zeta(npt))
      zeta = POINT_ZETA
      ! Unscreened (threshold 0): a one-shot analysis, where a pair the
      ! SCF's screen would drop costs nothing worth saving.
      call pairs%build(b, 0.0_dp)
      call pairs%to_device()
      !$acc enter data copyin(pts, zeta, density) create(electronic)
      call trc_pcm_potential(b, pairs, npt, pts, zeta, density, electronic)
      !$acc update self(electronic)
      !$acc exit data delete(pts, zeta, density, electronic)
      call pairs%release()
   end subroutine electronic_potential

   pure subroutine lattice_indices(g, nx1, nx2, i, j, k)
      !! Lattice point `g` (1-based, x fastest) as 0-based (i, j, k)
      !$acc routine seq
      integer, intent(in) :: g, nx1, nx2
      integer, intent(out) :: i, j, k
      integer :: rest
      rest = g - 1
      i = mod(rest, nx1)
      rest = rest/nx1
      j = mod(rest, nx2)
      k = rest/nx2
   end subroutine lattice_indices

   !
   ! KERNELS. Each is one `do concurrent` over a single flattened index,
   ! calling a same-module `pure` `!$acc routine seq` item routine with
   ! explicit-shape dummies and its scratch declared inside. Each item's sum
   ! over grid points runs in the points' own order, so a result does not
   ! depend on how the items are scheduled (host and device agree to the
   ! FMA contraction).
   !

   subroutine chelpg_lattice_kernel(n_lattice, nx1, nx2, low, h, margin, natm, at_r, radii, kept)
      !! Whether each lattice point survives the two exclusion tests
      integer, intent(in) :: n_lattice, nx1, nx2, natm
      real(dp), intent(in) :: low(3)
      real(dp), intent(in) :: h, margin
      real(dp), intent(in) :: at_r(3, natm), radii(natm)
      logical, intent(out) :: kept(n_lattice)
      integer :: g

      !$acc enter data copyin(low, at_r, radii) create(kept)
      do concurrent(g=1:n_lattice)
         kept(g) = lattice_point_kept(g, nx1, nx2, low, h, margin, natm, at_r, radii)
      end do
      !$acc update self(kept)
      !$acc exit data delete(low, at_r, radii, kept)
   end subroutine chelpg_lattice_kernel

   pure function lattice_point_kept(g, nx1, nx2, low, h, margin, natm, at_r, radii) result(keep)
      !$acc routine seq
      integer, intent(in) :: g, nx1, nx2, natm
      real(dp), intent(in) :: low(3)
      real(dp), intent(in) :: h, margin
      real(dp), intent(in) :: at_r(3, natm), radii(natm)
      logical :: keep

      real(dp) :: r(3)
      real(dp) :: dx, dy, dz, dist, nearest
      integer :: i, j, k, ia

      call lattice_indices(g, nx1, nx2, i, j, k)
      r(1) = low(1) + i*h
      r(2) = low(2) + j*h
      r(3) = low(3) + k*h
      keep = .false.
      nearest = huge(1.0_dp)
      do ia = 1, natm
         dx = r(1) - at_r(1, ia)
         dy = r(2) - at_r(2, ia)
         dz = r(3) - at_r(3, ia)
         dist = sqrt(dx*dx + dy*dy + dz*dz)
         if (dist < radii(ia)) return
         nearest = min(nearest, dist)
      end do
      keep = nearest <= margin
   end function lattice_point_kept

   subroutine chelpg_nuclear_kernel(natm, at_r, at_z, npt, pts, electronic, potential)
      !! potential(g) = sum_A Z_A/|r_g - R_A| - electronic(g)
      integer, intent(in) :: natm, npt
      real(dp), intent(in) :: at_r(3, natm), at_z(natm)
      real(dp), intent(in) :: pts(3, npt), electronic(npt)
      real(dp), intent(out) :: potential(npt)
      integer :: g

      !$acc enter data copyin(at_r, at_z, pts, electronic) create(potential)
      do concurrent(g=1:npt)
         potential(g) = nuclear_point(g, natm, at_r, at_z, npt, pts) - electronic(g)
      end do
      !$acc update self(potential)
      !$acc exit data delete(at_r, at_z, pts, electronic, potential)
   end subroutine chelpg_nuclear_kernel

   pure function nuclear_point(g, natm, at_r, weights, npt, pts) result(v)
      !! sum_A weights_A/|r_g - R_A|: the nuclei's potential, or the fitted charges'
      !$acc routine seq
      integer, intent(in) :: g, natm, npt
      real(dp), intent(in) :: at_r(3, natm), weights(natm), pts(3, npt)
      real(dp) :: v
      real(dp) :: dx, dy, dz
      integer :: ia
      v = 0.0_dp
      do ia = 1, natm
         dx = pts(1, g) - at_r(1, ia)
         dy = pts(2, g) - at_r(2, ia)
         dz = pts(3, g) - at_r(3, ia)
         v = v + weights(ia)/sqrt(dx*dx + dy*dy + dz*dz)
      end do
   end function nuclear_point

   subroutine chelpg_normal_kernel(natm, at_r, npt, pts, potential, gram, projection)
      !! gram(A + natm (B-1)) = sum_g 1/(r_gA r_gB) and projection(A) = sum_g V(g)/r_gA
      integer, intent(in) :: natm, npt
      real(dp), intent(in) :: at_r(3, natm), pts(3, npt), potential(npt)
      real(dp), intent(out) :: gram(natm*natm), projection(natm)
      integer :: ab

      ! One item per (A, B) and per A; O(natm^2 npt) work, the fit's
      ! dominant cost, and wide enough to fill the device from a few atoms.
      !$acc enter data copyin(at_r, pts, potential) create(gram, projection)
      do concurrent(ab=1:natm*natm)
         gram(ab) = gram_element(ab, natm, at_r, npt, pts)
      end do
      do concurrent(ab=1:natm)
         projection(ab) = projection_element(ab, natm, at_r, npt, pts, potential)
      end do
      !$acc update self(gram, projection)
      !$acc exit data delete(at_r, pts, potential, gram, projection)
   end subroutine chelpg_normal_kernel

   pure function gram_element(ab, natm, at_r, npt, pts) result(s)
      !$acc routine seq
      integer, intent(in) :: ab, natm, npt
      real(dp), intent(in) :: at_r(3, natm), pts(3, npt)
      real(dp) :: s
      real(dp) :: ax, ay, az, bx, by, bz
      integer :: g, ia, ib
      ia = 1 + mod(ab - 1, natm)
      ib = 1 + (ab - 1)/natm
      s = 0.0_dp
      do g = 1, npt
         ax = pts(1, g) - at_r(1, ia)
         ay = pts(2, g) - at_r(2, ia)
         az = pts(3, g) - at_r(3, ia)
         bx = pts(1, g) - at_r(1, ib)
         by = pts(2, g) - at_r(2, ib)
         bz = pts(3, g) - at_r(3, ib)
         s = s + (1.0_dp/sqrt(ax*ax + ay*ay + az*az))*(1.0_dp/sqrt(bx*bx + by*by + bz*bz))
      end do
   end function gram_element

   pure function projection_element(ia, natm, at_r, npt, pts, potential) result(s)
      !$acc routine seq
      integer, intent(in) :: ia, natm, npt
      real(dp), intent(in) :: at_r(3, natm), pts(3, npt), potential(npt)
      real(dp) :: s
      real(dp) :: ax, ay, az
      integer :: g
      s = 0.0_dp
      do g = 1, npt
         ax = pts(1, g) - at_r(1, ia)
         ay = pts(2, g) - at_r(2, ia)
         az = pts(3, g) - at_r(3, ia)
         s = s + potential(g)/sqrt(ax*ax + ay*ay + az*az)
      end do
   end function projection_element

   subroutine chelpg_residual_kernel(natm, at_r, charges, npt, pts, potential, residual)
      !! residual(g) = V(g) - sum_A q_A/|r_g - R_A|
      integer, intent(in) :: natm, npt
      real(dp), intent(in) :: at_r(3, natm), charges(natm), pts(3, npt), potential(npt)
      real(dp), intent(out) :: residual(npt)
      integer :: g

      !$acc enter data copyin(at_r, charges, pts, potential) create(residual)
      do concurrent(g=1:npt)
         residual(g) = potential(g) - nuclear_point(g, natm, at_r, charges, npt, pts)
      end do
      !$acc update self(residual)
      !$acc exit data delete(at_r, charges, pts, potential, residual)
   end subroutine chelpg_residual_kernel

end module trc_chelpg
