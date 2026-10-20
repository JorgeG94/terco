module trc_dft_partition
   !! Becke partition of space into atomic cells
   !!
   !! A molecular integral is done as a sum of atom-centred quadratures, which
   !! requires a partition of unity: a set of weights w_A(r) summing to 1 at
   !! every point, each concentrated near its atom. Becke's construction
   !! [JCP 88, 2547 (1988)] builds them from the confocal elliptical coordinate
   !! of each atom pair,
   !!
   !!    mu_AB = (|r - R_A| - |r - R_B|) / R_AB
   !!
   !! smoothed by a cutoff s(mu) that runs from 1 to 0 as mu goes -1 to 1. The
   !! unnormalised cell function is the product over partners,
   !!
   !!    P_A(r) = prod_{B /= A} s(nu_AB),   w_A = P_A / sum_C P_C
   !!
   !! Two cutoffs are offered. Becke's is three iterations of p(x) = x(3-x^2)/2,
   !! which is smooth everywhere and non-zero everywhere -- so every atom
   !! contributes at every point, and the cost is unavoidably natoms^2 per
   !! point. Stratmann's [CPL 257, 213 (1996)] is a degree-5 polynomial that
   !! reaches exactly +-1 at |mu| >= a, so distant pairs contribute exactly 0 or
   !! 1 and can be skipped. That is what makes large systems affordable.
   !!
   !! `nu` is `mu` after an adjustment for the two atoms having different sizes,
   !! so that the cell boundary sits between them rather than halfway. Becke's
   !! appendix gives the shift in terms of the radius ratio chi:
   !!
   !!    u = (chi - 1)/(chi + 1),   a = u/(u^2 - 1),  clamped to |a| <= 1/2
   !!
   !! which reduces to a = (1/chi - chi)/4. Treutler's variant uses the square
   !! roots of the radii instead. Both are available; the choice changes the
   !! weights slightly and must match whatever a reference is being compared to.
   use pic_types, only: dp, int64
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_dft_radial, only: bragg_radius
   implicit none
   private

   public :: PARTITION_BECKE, PARTITION_STRATMANN
   public :: ADJUST_NONE, ADJUST_BECKE, ADJUST_TREUTLER
   public :: becke_partition_weights
   public :: becke_partition_derivatives
   public :: becke_partition_grad_contract
   public :: partition_scheme_name
   public :: PARTITION_NU_MAX

   real(dp), parameter :: PARTITION_NU_MAX = 0.95_dp
   !! Screening threshold on the cutoff argument against the owner; see
   !! becke_partition_weights. Becke's s(0.95) is 1e-9; a cell function is a
   !! product of such factors, so a dropped atom's cell is at most that.

   ! Cutoff profile
   integer, parameter :: PARTITION_BECKE = 1     !! Three iterations of p(x), smooth everywhere
   integer, parameter :: PARTITION_STRATMANN = 2  !! Degree-5 with a hard cutoff, screenable

   ! Atomic size adjustment
   integer, parameter :: ADJUST_NONE = 0      !! All atoms the same size
   integer, parameter :: ADJUST_BECKE = 1     !! Bragg radii
   integer, parameter :: ADJUST_TREUTLER = 2  !! Square roots of the Bragg radii

   real(dp), parameter :: STRATMANN_A = 0.64_dp
   !! Stratmann's cutoff parameter, eq. 14

   real(dp), parameter :: MAX_ADJUST = 0.5_dp
   !! Becke's clamp on the size-adjustment shift

   real(dp), parameter, public :: TINY_RADIUS = 1.0e-200_dp
   !! Guards the ratio when a radius is zero (a ghost atom)

   real(dp), parameter :: CELL_FLOOR = 1.0e-150_dp
   !! A cell function below this is dropped from the gradient's sums in
   !! becke_partition_grad_contract (see partition_grad_point): every term
   !! it carries is below it, relatively, and above it every factor is
   !! safely non-zero to divide by.
   real(dp), parameter :: CELL_REL = 1.0e-20_dp
   !! ... and so is one below this times the point owner's own cell (see
   !! partition_grad_point for the bound on what that drops).

   integer, parameter :: N_DIM = 3

contains

   pure function partition_scheme_name(scheme) result(name)
      !! Human-readable scheme name, for logs and error messages
      integer, intent(in) :: scheme
      character(len=:), allocatable :: name

      select case (scheme)
      case (PARTITION_BECKE)
         name = "Becke"
      case (PARTITION_STRATMANN)
         name = "Stratmann"
      case default
         name = "unknown"
      end select
   end function partition_scheme_name

   subroutine becke_partition_weights(points, atom_coords, atomic_numbers, owner, &
                                      scheme, adjust, weights, error, nu_max, budget, first, last)
      !! Partition weight at each grid point, for the atom that owns it
      !!
      !! The returned weight is w_owner(r): the fraction of the integrand at
      !! that point assigned to the atom whose atomic grid produced it. Multiply
      !! it by the point's radial and angular weight to get its contribution.
      !!
      !! ON THE DEVICE. Every point is independent and the work per point is
      !! O(active atoms squared), which on a 74-atom molecule was seventeen
      !! seconds of one host thread before the first kernel ran. The loop is
      !! a `do concurrent` over points with the body in a `routine seq`
      !! procedure, terco's pattern. What made that awkward was the per-point
      !! scratch sized by the atom count -- distances, cell functions, the
      !! active list -- which a device thread cannot allocate. It lives in
      !! device arrays of (n_atoms, chunk) instead, and the points are
      !! processed chunk by chunk under `budget` doubles of scratch, so the
      !! atom count is bounded by memory and not by a parameter.
      real(dp), intent(in) :: points(:, :)       !! (3, n_points), Bohr
      real(dp), intent(in) :: atom_coords(:, :)  !! (3, n_atoms), Bohr
      integer, intent(in) :: atomic_numbers(:)   !! Z per atom, for the radii
      integer, intent(in) :: owner(:)            !! Atom each point belongs to
      integer, intent(in) :: scheme              !! PARTITION_BECKE or PARTITION_STRATMANN
      integer, intent(in) :: adjust              !! ADJUST_NONE, _BECKE or _TREUTLER
      real(dp), intent(out) :: weights(:)        !! (n_points)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: nu_max
      !! Screening. An atom whose cutoff argument against the point's owner
      !! is already at or beyond this value is left out of the pair loop for
      !! that point: with Becke's cutoff its cell function is then below
      !! 1e-9 relative, with Stratmann's it is exactly zero past 0.64. The
      !! loop is O(active^2) per point instead of O(n_atoms^2). Pass a value
      !! >= 1 for no screening; check_grid measures the difference. Default
      !! PARTITION_NU_MAX, or STRATMANN_A under that scheme.
      integer(int64), intent(in), optional :: budget
      !! Scratch budget in doubles, default 2^26 (512 MB): sets how many
      !! points a chunk holds, 2^26 / (4 n_atoms).
      integer, intent(in), optional :: first, last
      !! Only the points first..last are partitioned; every other weight is
      !! returned as zero, so a caller holding one block per rank can sum
      !! the arrays and have the whole grid. Default all of them.

      real(dp), allocatable :: shift(:, :), inv_distance(:, :), r_min(:)
      real(dp), allocatable :: atom_r(:, :), ract(:, :), cell(:, :)
      integer, allocatable :: active(:, :), mark(:, :)
      integer :: n_atoms, n_points, i, j, k0, n, chunk, g, k_lo, k_hi
      integer(int64) :: cap
      real(dp) :: numax

      n_atoms = size(atom_coords, 2)
      n_points = size(points, 2)

      if (size(atomic_numbers) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "partition: atomic_numbers does not match atom_coords")
         return
      end if
      if (size(owner) /= n_points .or. size(weights) /= n_points) then
         call error%set(ERROR_VALIDATION, "partition: owner and weights must match the points")
         return
      end if
      if (scheme /= PARTITION_BECKE .and. scheme /= PARTITION_STRATMANN) then
         call error%set(ERROR_VALIDATION, "partition: unknown scheme")
         return
      end if
      if (any(owner < 1) .or. any(owner > n_atoms)) then
         call error%set(ERROR_VALIDATION, "partition: owner index outside the atom list")
         return
      end if

      k_lo = 1
      k_hi = n_points
      if (present(first)) k_lo = max(first, 1)
      if (present(last)) k_hi = min(last, n_points)
      weights = 0.0_dp

      ! A single atom owns everything, and the pair loop below would leave the
      ! product empty. Short-circuit rather than special-case inside the loop.
      if (n_atoms == 1) then
         weights(k_lo:k_hi) = 1.0_dp
         return
      end if

      numax = PARTITION_NU_MAX
      ! Stratmann's cutoff is exactly zero from STRATMANN_A on, so there the
      ! screening is exact and nothing tighter buys anything.
      if (scheme == PARTITION_STRATMANN) numax = STRATMANN_A
      if (present(nu_max)) numax = nu_max
      cap = 67108864_int64
      if (present(budget)) cap = max(budget, int(4*n_atoms, int64))
      chunk = int(min(int(n_points, int64), max(1_int64, cap/int(4*n_atoms, int64))))

      allocate (shift(n_atoms, n_atoms), inv_distance(n_atoms, n_atoms), r_min(n_atoms))
      call size_adjustment(atomic_numbers, adjust, shift)

      ! 1/R_AB, precomputed: it is the same for every point.
      inv_distance = 0.0_dp
      do i = 1, n_atoms
         do j = 1, n_atoms
            if (i /= j) then
               inv_distance(i, j) = 1.0_dp/norm2(atom_coords(:, i) - atom_coords(:, j))
            end if
         end do
      end do

      ! Nearest neighbour of each atom, for Stratmann's inner sphere: a point
      ! closer to its owner than (1 - a)/2 of that distance has mu <= -a
      ! against every other atom, and the weight is exactly one.
      r_min = huge(1.0_dp)
      do i = 1, n_atoms
         do j = 1, n_atoms
            if (i /= j) r_min(i) = min(r_min(i), 1.0_dp/inv_distance(i, j))
         end do
      end do

      ! Point index FIRST: consecutive threads are consecutive points, and
      ! this is what makes their reads of the same column adjacent. `atom_r`
      ! and `mark` are indexed by atom, read with the atom loop uniform
      ! across threads; `ract` and `cell` are indexed by POSITION in each
      ! point's active list, so the O(active^2) pair loop -- whose atom
      ! indices differ from thread to thread -- still reads and writes
      ! whole columns. Indexing those two by atom made every pair two
      ! scattered read-modify-writes, and the partition memory-bound.
      allocate (atom_r(chunk, n_atoms), ract(chunk, n_atoms), cell(chunk, n_atoms), &
                active(chunk, n_atoms), mark(chunk, n_atoms))
      mark = 0
      ! Only the block this call writes is on the device: a copyout of the
      ! whole array would bring back uninitialised device memory for the
      ! rest, over the zeros set above.
      !$acc data copyin(points, atom_coords, owner, shift, inv_distance, r_min, mark) &
      !$acc      create(atom_r, ract, cell, active) copyout(weights(k_lo:k_hi))
      do k0 = k_lo, k_hi, chunk
         n = min(chunk, k_hi - k0 + 1)
         do concurrent(g=1:n)
            call partition_point(k0 + g - 1, g, n_atoms, n_points, chunk, points, atom_coords, owner, &
                                 shift, inv_distance, r_min, scheme, adjust, numax, &
                                 atom_r, ract, cell, active, mark, weights)
         end do
      end do
      !$acc end data
   end subroutine becke_partition_weights

   !
   ! One point's partition weight, into weights(k), using row g of the
   ! scratch. The atoms this point can see are first those not cut off
   ! against its owner, then any atom not cut off against one of THOSE,
   ! because an atom with a negligible cell of its own can still carry a
   ! factor of order one on a neighbour's cell (water: the far hydrogen is
   ! cut off against the oxygen and multiplies the near hydrogen's cell by
   ! 0.97). An atom cut off against every kept atom has a cell below
   ! s(nu_max) and a factor within s(nu_max) of one on all of them. `mark`
   ! is a membership mask, left all-zero on exit for the next point.
   !
   pure subroutine partition_point(k, g, n_atoms, n_points, chunk, points, atom_coords, owner, &
                                   shift, inv_distance, r_min, scheme, adjust, numax, &
                                   atom_r, ract, cell, active, mark, weights)
      !$acc routine seq
      integer, intent(in) :: k, g, n_atoms, n_points, chunk, scheme, adjust
      real(dp), intent(in) :: points(3, n_points), atom_coords(3, n_atoms)
      integer, intent(in) :: owner(n_points)
      real(dp), intent(in) :: shift(n_atoms, n_atoms), inv_distance(n_atoms, n_atoms), r_min(n_atoms)
      real(dp), intent(in) :: numax
      real(dp), intent(inout) :: atom_r(chunk, n_atoms), ract(chunk, n_atoms), cell(chunk, n_atoms)
      integer, intent(inout) :: active(chunk, n_atoms), mark(chunk, n_atoms)
      real(dp), intent(inout) :: weights(n_points)
      integer :: own, i, j, ia, ja, n_active, n1
      real(dp) :: mu, nu, s, total, r_inner, dx, dy, dz, ri, ci

      own = owner(k)
      do i = 1, n_atoms
         dx = points(1, k) - atom_coords(1, i)
         dy = points(2, k) - atom_coords(2, i)
         dz = points(3, k) - atom_coords(3, i)
         atom_r(g, i) = sqrt(dx*dx + dy*dy + dz*dz)
      end do
      if (scheme == PARTITION_STRATMANN .and. adjust == ADJUST_NONE) then
         r_inner = 0.5_dp*(1.0_dp - STRATMANN_A)*r_min(own)
         if (atom_r(g, own) < r_inner) then
            weights(k) = 1.0_dp
            return
         end if
      end if

      n_active = 1
      active(g, 1) = own
      ract(g, 1) = atom_r(g, own)
      mark(g, own) = 1
      do i = 1, n_atoms
         if (i == own) cycle
         mu = (atom_r(g, i) - atom_r(g, own))*inv_distance(i, own)
         nu = mu + shift(i, own)*(1.0_dp - mu*mu)
         if (nu < numax) then
            n_active = n_active + 1
            active(g, n_active) = i
            ract(g, n_active) = atom_r(g, i)
            mark(g, i) = 1
         end if
      end do
      n1 = n_active
      do i = 1, n_atoms
         if (mark(g, i) == 1) cycle
         ri = atom_r(g, i)
         do ja = 2, n1
            j = active(g, ja)
            mu = (ri - ract(g, ja))*inv_distance(i, j)
            nu = mu + shift(i, j)*(1.0_dp - mu*mu)
            if (nu < numax) then
               n_active = n_active + 1
               active(g, n_active) = i
               ract(g, n_active) = ri
               mark(g, i) = 1
               exit
            end if
         end do
      end do

      ! `cell` by active position: cell(g, ia) is atom active(g, ia)'s cell.
      ! Atom i's cell takes its (1 - s) factors from the earlier rows (as
      ! somebody's j) before its own row multiplies in its s factors, the
      ! same order as indexing by atom did, so the weights are unchanged
      ! to the last bit; the row's own product just runs in a register.
      do ia = 1, n_active
         cell(g, ia) = 1.0_dp
      end do
      do ia = 1, n_active
         i = active(g, ia)
         ri = ract(g, ia)
         ci = cell(g, ia)
         do ja = ia + 1, n_active
            j = active(g, ja)
            mu = (ri - ract(g, ja))*inv_distance(i, j)
            nu = mu + shift(i, j)*(1.0_dp - mu*mu)
            if (scheme == PARTITION_BECKE) then
               s = becke_cutoff(nu)
            else
               s = stratmann_cutoff(nu)
            end if
            ! s is the fraction going to i; the complement goes to j. Doing
            ! both here keeps the pair symmetric by construction.
            ci = ci*s
            cell(g, ja) = cell(g, ja)*(1.0_dp - s)
         end do
         cell(g, ia) = ci
      end do

      total = 0.0_dp
      do ia = 1, n_active
         total = total + cell(g, ia)
         mark(g, active(g, ia)) = 0
      end do
      if (total > 0.0_dp) then
         weights(k) = cell(g, 1)/total
      else
         ! Every cell function underflowed, which happens only absurdly far
         ! from the molecule where the integrand is zero anyway.
         weights(k) = 0.0_dp
      end if
   end subroutine partition_point

   subroutine becke_partition_derivatives(points, atom_coords, atomic_numbers, owner, &
                                          scheme, adjust, dweights, error, nu_max)
      !! dP_k/dR_A at fixed grid points, as (3, n_atoms, n_points)
      !!
      !! The partition weight depends on every nuclear position, so moving any
      !! atom changes the weight of every point -- including points owned by
      !! other atoms. That is the whole content of the "grid response" term in
      !! a DFT gradient, and omitting it leaves a gradient that looks entirely
      !! reasonable and disagrees with finite differences at about 1e-4.
      !!
      !! **The grid points are held fixed here.** Points move with the atom that
      !! owns them, but that motion is a separate contribution which the caller
      !! adds; mixing the two into one routine would make neither checkable.
      !!
      !! SCREENED THE SAME WAY `becke_partition_weights` IS. The energy this
      !! feeds does not sum the UNSCREENED partition product over every
      !! atom: `apply_partition` calls `becke_partition_weights` at its
      !! default `nu_max` (PARTITION_NU_MAX, or STRATMANN_A under that
      !! scheme), which drops an atom from a point's active set once its
      !! cutoff against that point is already below ~1e-9 relative, exactly
      !! `becke_partition_weights`'s own header. A gradient of the UNSCREENED
      !! partition is not the exact derivative of the SCREENED energy that
      !! actually ran -- a difference `partition_point`'s own screening
      !! keeps below the energy's own truncation, so it never mattered for
      !! the FD/pyscf numbers trc_xc_grad.F90's header reports on water-
      !! and NH2-sized cases, but it is the difference between O(active^2)
      !! and O(n_atoms^2) per point on anything larger, and this routine had
      !! no screening at all before. The active-set build below (own atom,
      !! then direct, then indirect cutoffs) is `partition_point`'s own,
      !! unchanged, so the SAME atoms the energy keeps are the only ones
      !! this derivative ever touches -- an atom the energy drops contributes
      !! exactly zero here too, which is correct under the screening
      !! approximation (the screened weight literally does not depend on
      !! that atom's position), not an extra error on top of it.
      !!
      !! O(active_atoms^2) PER POINT, O(n_atoms) scratch per point -- an
      !! EARLIER version of this routine was unscreened O(n_atoms^3) (see
      !! git history and check_dft_partition_deriv.F90's reference copy).
      !!
      !! cell_i = prod_{j /= i, j active} S(i,j), S(i,j) the fraction of the
      !! (i,j) pair that goes to i (S(i,j) + S(j,i) = 1, from
      !! `becke_cutoff`/`stratmann_cutoff` of the size-adjusted confocal
      !! coordinate). Every factor S(i,j) depends on R_i and R_j ONLY -- so
      !! d(cell_i)/dR_a for a /= i picks out exactly the ONE factor touching
      !! a:
      !!
      !!     d(cell_i)/dR_a = [dS(i,a)/dR_a] * prod_{k /= i,a, active} S(i,k)     (a /= i)
      !!     d(cell_i)/dR_i = sum_{j /= i, active} [dS(i,j)/dR_i] * prod_{k /= i,j, active} S(i,k)
      !!
      !! "prod_{k/=i,a} S(i,k)" is cell_i with ONE factor excluded, for every
      !! active a in turn -- the classic prefix/suffix-product trick (never a
      !! division, so a cutoff factor underflowing to exactly zero, which
      !! value and derivative both do legitimately, is never a 0/0): for row
      !! i, lay the active factors S(i,j) out in a list, take running
      !! products from the front (`prefix`) and from the back (`suffix`),
      !! and "all but the t-th" is `prefix(t-1)*suffix(t+1)`. One pass per
      !! active row i gives every active a /= i AND the i = a term in the
      !! same sweep, so the active-by-active derivative matrix this routine
      !! ultimately needs (every d(cell_i)/dR_a, because `dtotal` below sums
      !! over every active i) costs O(active) work `n_active` times --
      !! O(active^2) total, with nothing sized n_atoms^2 kept around:
      !! `dtotal`/`d_own` are O(n_atoms) (indexed by real atom, zero outside
      !! the active set) and the S(i,j)/dS(i,j) values are recomputed on the
      !! fly from `atom_r`/`shift`/`inv_distance` rather than stored.
      !!
      !! Checked against an UNSCREENED reference implementation on water
      !! and formamide -- see check_dft_partition_deriv.F90 -- to ~1e-14 (at
      !! `nu_max >= 1`, screening off on both sides; the two are NOT expected
      !! to agree at the production `nu_max`, by design), and the active-set
      !! build itself against `becke_partition_weights`'s FD (perturbing one
      !! atom and differencing the SCREENED weight) on a case large enough
      !! for screening to actually drop atoms -- see
      !! check_xc_grad_fd.F90's own larger case.
      real(dp), intent(in) :: points(:, :)       !! (3, n_points), Bohr
      real(dp), intent(in) :: atom_coords(:, :)  !! (3, n_atoms), Bohr
      integer, intent(in) :: atomic_numbers(:)   !! Z per atom, for the radii
      integer, intent(in) :: owner(:)            !! Atom each point belongs to
      integer, intent(in) :: scheme
      integer, intent(in) :: adjust
      real(dp), intent(out) :: dweights(:, :, :)  !! (3, n_atoms, n_points)
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: nu_max
      !! As `becke_partition_weights`'s own `nu_max`: default PARTITION_NU_MAX,
      !! or STRATMANN_A under that scheme -- the SAME default the energy's own
      !! unqualified `becke_partition_weights` call (`apply_partition`) uses.
      !! Pass >= 1 for no screening (what check_dft_partition_deriv.F90 does
      !! to compare against the unscreened reference).

      real(dp), allocatable :: shift(:, :), inv_distance(:, :), r_min(:)
      real(dp), allocatable :: atom_r(:), unit_vec(:, :)
      real(dp), allocatable :: f(:), prefix(:), suffix(:)
      real(dp), allocatable :: own_deriv(:, :), oth_deriv(:, :)
      real(dp), allocatable :: dtotal(:, :), d_own(:, :)
      integer, allocatable :: idx(:), active(:), mark(:)
      real(dp) :: cell_own, total, excl, numax, r_inner, mu, nu
      real(dp) :: d_ii(3), d_ia(3)
      integer :: n_atoms, n_points, i, j, k, t, m, own, ia, ja, n_active, n1

      n_atoms = size(atom_coords, 2)
      n_points = size(points, 2)
      dweights = 0.0_dp

      if (size(atomic_numbers) /= n_atoms) then
         call error%set(ERROR_VALIDATION, &
                        "partition derivatives: atomic_numbers does not match atom_coords")
         return
      end if
      if (scheme /= PARTITION_BECKE .and. scheme /= PARTITION_STRATMANN) then
         call error%set(ERROR_VALIDATION, "partition derivatives: unknown scheme")
         return
      end if

      ! One atom owns everything and its weight is identically one, so it does
      ! not change when the atom moves.
      if (n_atoms == 1) return

      numax = PARTITION_NU_MAX
      if (scheme == PARTITION_STRATMANN) numax = STRATMANN_A
      if (present(nu_max)) numax = nu_max

      allocate (shift(n_atoms, n_atoms), inv_distance(n_atoms, n_atoms), r_min(n_atoms))

      call size_adjustment(atomic_numbers, adjust, shift)

      inv_distance = 0.0_dp
      do i = 1, n_atoms
         do j = 1, n_atoms
            if (i /= j) then
               inv_distance(i, j) = 1.0_dp/norm2(atom_coords(:, i) - atom_coords(:, j))
            end if
         end do
      end do

      ! Nearest neighbour of each atom, for Stratmann's inner sphere -- the
      ! same shortcut `partition_point` takes (a point closer to its owner
      ! than that has weight exactly 1, so EVERY atom's derivative there is
      ! exactly zero, not merely small: `dweights` stays at the zero this
      ! routine already initialises it to, via the `cycle` below).
      r_min = huge(1.0_dp)
      do i = 1, n_atoms
         do j = 1, n_atoms
            if (i /= j) r_min(i) = min(r_min(i), 1.0_dp/inv_distance(i, j))
         end do
      end do

      ! One thread per grid point, as before: every iteration reads the
      ! geometry and writes its own `dweights(:, :, k)`, so the points are
      ! independent. `shift`/`inv_distance`/`r_min` depend only on the
      ! nuclei and stay shared; everything else is O(n_atoms) scratch,
      ! allocated once per thread.
      !$omp parallel default(none) &
      !$omp    shared(points, atom_coords, owner, scheme, dweights, shift, &
      !$omp           inv_distance, r_min, numax, n_points, n_atoms, adjust) &
      !$omp    private(k, i, j, t, m, own, ia, ja, n_active, n1, cell_own, total, excl, &
      !$omp            r_inner, mu, nu, d_ii, d_ia, atom_r, unit_vec, f, prefix, suffix, &
      !$omp            own_deriv, oth_deriv, dtotal, d_own, idx, active, mark)
      allocate (atom_r(n_atoms), unit_vec(3, n_atoms))
      allocate (f(n_atoms), prefix(0:n_atoms), suffix(n_atoms + 1), idx(n_atoms))
      allocate (own_deriv(3, n_atoms), oth_deriv(3, n_atoms))
      allocate (dtotal(3, n_atoms), d_own(3, n_atoms))
      allocate (active(n_atoms), mark(n_atoms))

      !$omp do schedule(static)
      do k = 1, n_points
         own = owner(k)
         do i = 1, n_atoms
            atom_r(i) = norm2(points(:, k) - atom_coords(:, i))
            if (atom_r(i) > 0.0_dp) then
               unit_vec(:, i) = (points(:, k) - atom_coords(:, i))/atom_r(i)
            else
               ! The point sits on the nucleus. |r - R| is not differentiable
               ! there, and the quadrature never places a point exactly on one,
               ! so this is a guard rather than a case.
               unit_vec(:, i) = 0.0_dp
            end if
         end do

         if (scheme == PARTITION_STRATMANN .and. adjust == ADJUST_NONE) then
            r_inner = 0.5_dp*(1.0_dp - STRATMANN_A)*r_min(own)
            if (atom_r(own) < r_inner) cycle
         end if

         ! The active set, exactly as `partition_point` builds it: the
         ! owner always first, then every atom not cut off against it,
         ! then (second pass) every remaining atom not cut off against any
         ! atom ALREADY active -- an atom with a negligible cell of its own
         ! can still carry a factor of order one on a neighbour's.
         mark = 0
         n_active = 1
         active(1) = own
         mark(own) = 1
         do i = 1, n_atoms
            if (i == own) cycle
            mu = (atom_r(i) - atom_r(own))*inv_distance(i, own)
            nu = mu + shift(i, own)*(1.0_dp - mu*mu)
            if (nu < numax) then
               n_active = n_active + 1
               active(n_active) = i
               mark(i) = 1
            end if
         end do
         n1 = n_active
         do i = 1, n_atoms
            if (mark(i) == 1) cycle
            do ja = 2, n1
               j = active(ja)
               mu = (atom_r(i) - atom_r(j))*inv_distance(i, j)
               nu = mu + shift(i, j)*(1.0_dp - mu*mu)
               if (nu < numax) then
                  n_active = n_active + 1
                  active(n_active) = i
                  mark(i) = 1
                  exit
               end if
            end do
         end do

         dtotal = 0.0_dp
         cell_own = 1.0_dp
         total = 0.0_dp

         do ia = 1, n_active
            i = active(ia)
            ! Row i: every ACTIVE factor of cell_i, in some fixed order,
            ! plus its two one-sided derivatives (own_deriv: d S(i,j)/dR_i;
            ! oth_deriv: d S(i,j)/dR_j) -- pair_sij below handles both
            ! orderings of (i, j) off the SAME canonical-pair (p < q)
            ! evaluation this routine always used.
            m = 0
            do ja = 1, n_active
               j = active(ja)
               if (j == i) cycle
               m = m + 1
               idx(m) = j
               call pair_sij(i, j, n_atoms, atom_coords, atom_r, unit_vec, shift, &
                            inv_distance, scheme, f(m), own_deriv(:, m), oth_deriv(:, m))
            end do

            prefix(0) = 1.0_dp
            do t = 1, m
               prefix(t) = prefix(t - 1)*f(t)
            end do
            suffix(m + 1) = 1.0_dp
            do t = m, 1, -1
               suffix(t) = suffix(t + 1)*f(t)
            end do

            d_ii = 0.0_dp
            do t = 1, m
               excl = prefix(t - 1)*suffix(t + 1)
               d_ii = d_ii + own_deriv(:, t)*excl
               j = idx(t)
               d_ia = oth_deriv(:, t)*excl
               dtotal(:, j) = dtotal(:, j) + d_ia
               if (i == own) d_own(:, j) = d_ia
            end do
            dtotal(:, i) = dtotal(:, i) + d_ii
            ! prefix(m) is cell_i = prod over ACTIVE j/=i of S(i,j) -- this
            ! row's own full product, already paid for building
            ! `d_ii`/`d_ia`, so `total` (sum over active i of cell_i) costs
            ! nothing extra here: no second pass to recover it.
            total = total + prefix(m)
            if (i == own) then
               d_own(:, i) = d_ii
               cell_own = prefix(m)
            end if
         end do

         if (total <= 0.0_dp) cycle

         do ia = 1, n_active
            i = active(ia)
            dweights(:, i, k) = (d_own(:, i)*total - cell_own*dtotal(:, i))/(total*total)
         end do
      end do
      !$omp end do
      !$omp end parallel
   end subroutine becke_partition_derivatives

   subroutine becke_partition_grad_contract(n_points, points, atom_coords, atomic_numbers, owner, &
                                            coef, scheme, adjust, grad, error, nu_max, budget)
      !! sum_k coef_k dP_k/dR_A, plus the point-motion term, straight into
      !! a (3, n_atoms) gradient
      !!
      !! What a DFT gradient's grid response needs from the partition is
      !! not dP_k/dR_A itself but its contraction with one number per point
      !! (the quadrature weight times the energy density), together with
      !! the motion of the point with its owner: P_k depends on relative
      !! positions only, so grad_{r_k} P_k = - sum_A dP_k/dR_A, landing on
      !! owner(k). Both are accumulated here, per point, into `grad`:
      !!
      !!     grad(:, A)        += coef_k dP_k/dR_A          (every active A)
      !!     grad(:, owner(k)) -= coef_k sum_A dP_k/dR_A
      !!
      !! which is what `becke_partition_derivatives` followed by that
      !! double loop gives, without its dense (3, n_atoms, n_points) array:
      !! that is O(n_atoms) memory per point, the limit on a large molecule.
      !!
      !! ON THE DEVICE, as `becke_partition_weights`: one point per thread,
      !! the per-point scratch in (chunk, ...) device arrays, point index
      !! first. Everything the O(active^2) pair loops touch (distances,
      !! factors, suffix products, derivative sums) is indexed by POSITION
      !! in the point's active list, not by atom, so a warp reads and
      !! writes whole columns instead of scattered entries; the running
      !! prefix product of each row is a register. The active set,
      !! screening and arithmetic order per point are those of
      !! `becke_partition_derivatives`, so each point's dP/dR is the same
      !! number; only the order the points' contributions are summed in
      !! differs (atomics on the device).
      !!
      !! `points`, `owner` and `coef` may already be on the device (the
      !! XC gradient's are); `grad` is a host array, ADDED into.
      integer, intent(in) :: n_points
      real(dp), intent(in) :: points(3, n_points)   !! Bohr
      real(dp), intent(in) :: atom_coords(:, :)     !! (3, n_atoms), Bohr
      integer, intent(in) :: atomic_numbers(:)      !! Z per atom, for the radii
      integer, intent(in) :: owner(n_points)        !! Atom each point belongs to
      real(dp), intent(in) :: coef(n_points)
      !! The per-point factor; a point with coef exactly 0 is skipped.
      integer, intent(in) :: scheme, adjust
      real(dp), intent(inout) :: grad(:, :)         !! (3, n_atoms), added into
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: nu_max
      !! As `becke_partition_derivatives`'s.
      integer(int64), intent(in), optional :: budget
      !! Scratch budget in doubles, default 2^26 (512 MB): a chunk holds
      !! 2^26 / (4 n_atoms) points, 93k at 180 atoms.

      real(dp), allocatable :: shift(:, :), inv_distance(:, :), r_min(:)
      real(dp), allocatable :: ract(:, :), cell(:, :)
      integer, allocatable :: active(:, :), sig(:, :), mark(:, :)
      integer :: n_atoms, i, j, k0, n, chunk, g
      integer(int64) :: cap
      real(dp) :: numax

      n_atoms = size(atom_coords, 2)
      if (size(atomic_numbers) /= n_atoms .or. size(grad, 2) /= n_atoms .or. size(grad, 1) /= 3) then
         call error%set(ERROR_VALIDATION, "partition gradient: atomic_numbers/grad do not match atom_coords")
         return
      end if
      if (scheme /= PARTITION_BECKE .and. scheme /= PARTITION_STRATMANN) then
         call error%set(ERROR_VALIDATION, "partition gradient: unknown scheme")
         return
      end if
      ! One atom owns everything: its weight is identically one.
      if (n_atoms == 1 .or. n_points == 0) return

      numax = PARTITION_NU_MAX
      if (scheme == PARTITION_STRATMANN) numax = STRATMANN_A
      if (present(nu_max)) numax = nu_max
      cap = 67108864_int64
      if (present(budget)) cap = max(budget, int(4*n_atoms, int64))
      chunk = int(min(int(n_points, int64), max(1_int64, cap/int(4*n_atoms, int64))))

      allocate (shift(n_atoms, n_atoms), inv_distance(n_atoms, n_atoms), r_min(n_atoms))
      call size_adjustment(atomic_numbers, adjust, shift)
      inv_distance = 0.0_dp
      do i = 1, n_atoms
         do j = 1, n_atoms
            if (i /= j) inv_distance(i, j) = 1.0_dp/norm2(atom_coords(:, i) - atom_coords(:, j))
         end do
      end do
      r_min = huge(1.0_dp)
      do i = 1, n_atoms
         do j = 1, n_atoms
            if (i /= j) r_min(i) = min(r_min(i), 1.0_dp/inv_distance(i, j))
         end do
      end do

      allocate (ract(chunk, n_atoms), cell(chunk, n_atoms))
      allocate (active(chunk, n_atoms), sig(chunk, n_atoms), mark(chunk, n_atoms))
      mark = 0
      !$acc data copyin(atom_coords, shift, inv_distance, r_min, mark) copy(grad) &
      !$acc      create(ract, cell, active, sig)
      do k0 = 1, n_points, chunk
         n = min(chunk, n_points - k0 + 1)
#ifdef _OPENACC
         do concurrent(g=1:n)
            call partition_grad_point(k0 + g - 1, g, n_atoms, n_points, chunk, points, atom_coords, owner, &
                                      coef, shift, inv_distance, r_min, scheme, adjust, numax, &
                                      ract, cell, active, sig, mark, grad)
         end do
#else
         ! Serial on the host: the body's atomics are inert there, and two
         ! points can land on the same atom.
         do g = 1, n
            call partition_grad_point(k0 + g - 1, g, n_atoms, n_points, chunk, points, atom_coords, owner, &
                                      coef, shift, inv_distance, r_min, scheme, adjust, numax, &
                                      ract, cell, active, sig, mark, grad)
         end do
#endif
      end do
      !$acc end data
   end subroutine becke_partition_grad_contract

   !
   ! One point of `becke_partition_grad_contract`, row g of the scratch.
   !
   ! The active set is `partition_point`'s. With cell_i = prod_{j active,
   ! j /= i} S(i,j) and the point's weight P = cell_own / total, total the
   ! sum of the active cells,
   !
   !     dP/dR_c = (dcell_own/dR_c total - cell_own dtotal/dR_c) / total^2
   !     dtotal/dR_c = sum_{r /= c} S'(r,c) cell_r/S(r,c)               (A)
   !                 - sum_{r /= c} S'(r,c) cell_c/S(c,r)               (B)
   !
   ! where S'(r,c) = dS(r,c)/dR_c, the only factor of cell_r that moves
   ! with R_c (A), and (B) is cell_c's own derivative, using S(c,r) =
   ! 1 - S(r,c). dcell_own/dR_c is the r = own term of (A) for c /= own
   ! and (B) for c = own. Everything for nucleus c is a sum over r, so one
   ! pass per c accumulates it in registers and adds it to the gradient
   ! straight away: the per-point scratch is the active list, its
   ! distances and the cells -- three columns, read the same way by every
   ! thread of a warp -- and nothing of size active^2 or 3 x active is
   ! read back and forth. That traffic was what bound the kernel.
   !
   ! "cell_r / S(r,c)" is the product of cell_r's other factors only while
   ! cell_r is a normal number, and most active atoms' cells are far below
   ! any that matters: a row whose cell is below CELL_FLOOR, or below
   ! CELL_REL times the owner's cell, counts as zero (its product stops as
   ! soon as it is there, and it leaves `total`). Its terms are that cell
   ! times S'/S, at most 8/(1 - nu) |grad nu| for Becke's cutoff -- under
   ! 1e-14 relative even at nu = 1 - 1e-6 -- and no factor that is exactly
   ! zero is ever divided by. A nucleus whose cell is zero then needs only
   ! (A), over the rows that are not (`sig`), so the pass over nuclei costs
   ! active x significant rather than active^2. The pair algebra is
   ! `pair_sij`'s, operation for operation, on scalars.
   !
   pure subroutine partition_grad_point(k, g, n_atoms, n_points, chunk, points, atom_coords, owner, &
                                        coef, shift, inv_distance, r_min, scheme, adjust, numax, &
                                        ract, cell, active, sig, mark, grad)
      !$acc routine seq
      integer, intent(in) :: k, g, n_atoms, n_points, chunk, scheme, adjust
      real(dp), intent(in) :: points(3, n_points), atom_coords(3, n_atoms)
      integer, intent(in) :: owner(n_points)
      real(dp), intent(in) :: coef(n_points)
      real(dp), intent(in) :: shift(n_atoms, n_atoms), inv_distance(n_atoms, n_atoms), r_min(n_atoms)
      real(dp), intent(in) :: numax
      real(dp), intent(inout) :: ract(chunk, n_atoms), cell(chunk, n_atoms)
      integer, intent(inout) :: active(chunk, n_atoms), sig(chunk, n_atoms), mark(chunk, n_atoms)
      real(dp), intent(inout) :: grad(3, n_atoms)
      integer :: own, i, j, ia, ja, n_active, n1, pp, qq, ic, ir, it, nsig, nrow
      real(dp) :: c, mu, nu, ri, rown, s, total, cell_own, cell_c, cell_r, r_inner, inv_total, w_own, floor
      real(dp) :: d3(3)
      real(dp) :: ucx, ucy, ucz, rc, invd, r_p, r_q, sh, ds_dnu, c1, t3, f0, f1, f2, f3
      real(dp) :: s_rc, s_cr, rx, ry, rz, dx, dy, dz, ex, ey
      real(dp) :: tax, tay, taz, tbx, tby, tbz, ox, oy, oz
      real(dp) :: vx, vy, vz, vox, voy, voz

      c = coef(k)
      if (c == 0.0_dp) return
      own = owner(k)
      ! Distances are recomputed where they are needed rather than kept
      ! per atom: less scratch per point is more points per launch.
      d3(1) = points(1, k) - atom_coords(1, own)
      d3(2) = points(2, k) - atom_coords(2, own)
      d3(3) = points(3, k) - atom_coords(3, own)
      rown = norm2(d3)
      if (scheme == PARTITION_STRATMANN .and. adjust == ADJUST_NONE) then
         ! Weight exactly one inside the inner sphere: no derivative.
         r_inner = 0.5_dp*(1.0_dp - STRATMANN_A)*r_min(own)
         if (rown < r_inner) return
      end if

      n_active = 1
      active(g, 1) = own
      ract(g, 1) = rown
      mark(g, own) = 1
      do i = 1, n_atoms
         if (i == own) cycle
         d3(1) = points(1, k) - atom_coords(1, i)
         d3(2) = points(2, k) - atom_coords(2, i)
         d3(3) = points(3, k) - atom_coords(3, i)
         ri = norm2(d3)
         mu = (ri - rown)*inv_distance(i, own)
         nu = mu + shift(i, own)*(1.0_dp - mu*mu)
         if (nu < numax) then
            n_active = n_active + 1
            active(g, n_active) = i
            ract(g, n_active) = ri
            mark(g, i) = 1
         end if
      end do
      n1 = n_active
      do i = 1, n_atoms
         if (mark(g, i) == 1) cycle
         d3(1) = points(1, k) - atom_coords(1, i)
         d3(2) = points(2, k) - atom_coords(2, i)
         d3(3) = points(3, k) - atom_coords(3, i)
         ri = norm2(d3)
         do ja = 2, n1
            j = active(g, ja)
            mu = (ri - ract(g, ja))*inv_distance(i, j)
            nu = mu + shift(i, j)*(1.0_dp - mu*mu)
            if (nu < numax) then
               n_active = n_active + 1
               active(g, n_active) = i
               ract(g, n_active) = ri
               mark(g, i) = 1
               exit
            end if
         end do
      end do

      ! The cells, each in a register, in becke_partition_derivatives'
      ! factor order (j by position, skipping i). The owner's comes first
      ! (it is position 1); every other row stops as soon as its product,
      ! which only falls, is below `floor`, and counts as zero from then on.
      ! The rows that do not are the `sig` list.
      total = 0.0_dp
      cell_own = 0.0_dp
      floor = CELL_FLOOR
      nsig = 0
      do ia = 1, n_active
         i = active(g, ia)
         ri = ract(g, ia)
         cell_c = 1.0_dp
         do ja = 1, n_active
            if (ja == ia) cycle
            j = active(g, ja)
            pp = min(i, j)
            qq = max(i, j)
            if (i < j) then
               r_p = ri
               r_q = ract(g, ja)
            else
               r_p = ract(g, ja)
               r_q = ri
            end if
            mu = (r_p - r_q)*inv_distance(pp, qq)
            nu = mu + shift(pp, qq)*(1.0_dp - mu*mu)
            if (scheme == PARTITION_BECKE) then
               f0 = nu
               f1 = 0.5_dp*f0*(3.0_dp - f0*f0)
               f2 = 0.5_dp*f1*(3.0_dp - f1*f1)
               f3 = 0.5_dp*f2*(3.0_dp - f2*f2)
               s = 0.5_dp*(1.0_dp - f3)
            else
               s = stratmann_cutoff(nu)
            end if
            if (i < j) then
               cell_c = cell_c*s
            else
               cell_c = cell_c*(1.0_dp - s)
            end if
            if (cell_c < floor) exit
         end do
         if (ia == 1) then
            cell_own = cell_c
            floor = max(CELL_FLOOR, CELL_REL*cell_own)
         end if
         if (cell_c < floor) cell_c = 0.0_dp
         cell(g, ia) = cell_c
         if (cell_c > 0.0_dp) then
            nsig = nsig + 1
            sig(g, nsig) = ia
            total = total + cell_c
         end if
      end do

      if (total > 0.0_dp) then
         ! 1/total and P = cell_own/total, never total**2: a point whose
         ! cells all sit near CELL_FLOOR has total ~1e-150, total**2 is
         ! below the smallest normal double, and a flush-to-zero build
         ! (nvfortran -fast) then divides by an exact zero.
         inv_total = 1.0_dp/total
         w_own = cell_own*inv_total
         vox = 0.0_dp
         voy = 0.0_dp
         voz = 0.0_dp
         do ic = 1, n_active
            j = active(g, ic)
            rc = ract(g, ic)
            cell_c = cell(g, ic)
            if (rc > 0.0_dp) then
               ucx = (points(1, k) - atom_coords(1, j))/rc
               ucy = (points(2, k) - atom_coords(2, j))/rc
               ucz = (points(3, k) - atom_coords(3, j))/rc
            else
               ucx = 0.0_dp
               ucy = 0.0_dp
               ucz = 0.0_dp
            end if
            ! (A), (B), and dcell_own/dR_c, summed over the rows r: every
            ! row when c's own cell counts ((B) needs them all), only the
            ! `sig` rows when it does not ((A) is all that is left).
            tax = 0.0_dp
            tay = 0.0_dp
            taz = 0.0_dp
            tbx = 0.0_dp
            tby = 0.0_dp
            tbz = 0.0_dp
            ox = 0.0_dp
            oy = 0.0_dp
            oz = 0.0_dp
            if (cell_c > 0.0_dp) then
               nrow = n_active
            else
               nrow = nsig
            end if
            do it = 1, nrow
               if (cell_c > 0.0_dp) then
                  ir = it
               else
                  ir = sig(g, it)
               end if
               if (ir == ic) cycle
               cell_r = cell(g, ir)
               i = active(g, ir)
               ri = ract(g, ir)
               pp = min(i, j)
               qq = max(i, j)
               if (i < j) then
                  r_p = ri
                  r_q = rc
               else
                  r_p = rc
                  r_q = ri
               end if
               invd = inv_distance(pp, qq)
               sh = shift(pp, qq)
               rx = atom_coords(1, pp) - atom_coords(1, qq)
               ry = atom_coords(2, pp) - atom_coords(2, qq)
               rz = atom_coords(3, pp) - atom_coords(3, qq)
               mu = (r_p - r_q)*invd
               nu = mu + sh*(1.0_dp - mu*mu)
               if (scheme == PARTITION_BECKE) then
                  f0 = nu
                  f1 = 0.5_dp*f0*(3.0_dp - f0*f0)
                  f2 = 0.5_dp*f1*(3.0_dp - f1*f1)
                  f3 = 0.5_dp*f2*(3.0_dp - f2*f2)
                  s = 0.5_dp*(1.0_dp - f3)
                  ds_dnu = -0.5_dp*(1.5_dp*(1.0_dp - f2*f2)) &
                           *(1.5_dp*(1.0_dp - f1*f1)) &
                           *(1.5_dp*(1.0_dp - f0*f0))
               else
                  s = stratmann_cutoff(nu)
                  ds_dnu = stratmann_cutoff_derivative(nu)
               end if
               c1 = ds_dnu*(1.0_dp - 2.0_dp*sh*mu)
               t3 = (r_p - r_q)*invd**3
               ! S'(r,c) = dS(r,c)/dR_c, pair_sij's d_oth for (i, j) = (r, c):
               ! ds_q when r < c (c is q), -ds_p when r > c (c is p).
               if (i < j) then
                  s_rc = s
                  s_cr = 1.0_dp - s
                  dx = c1*(ucx*invd + t3*rx)
                  dy = c1*(ucy*invd + t3*ry)
                  dz = c1*(ucz*invd + t3*rz)
               else
                  s_rc = 1.0_dp - s
                  s_cr = s
                  dx = -(c1*(-ucx*invd - t3*rx))
                  dy = -(c1*(-ucy*invd - t3*ry))
                  dz = -(c1*(-ucz*invd - t3*rz))
               end if
               ! cell/S is the product of the cell's OTHER factors. Where S
               ! itself is exactly 0 (Becke's s rounds to 0 or 1 in double,
               ! and a -fast build contracts the polynomial into FMAs that
               ! reach it) the term is S'*(that product), and S' is zero
               ! there too (ds/dnu vanishes at nu = +-1): skip it rather
               ! than form 0/0.
               if (cell_r > 0.0_dp .and. s_rc > 0.0_dp) then
                  ex = cell_r/s_rc
                  tax = tax + dx*ex
                  tay = tay + dy*ex
                  taz = taz + dz*ex
                  if (ir == 1) then
                     ox = dx*ex
                     oy = dy*ex
                     oz = dz*ex
                  end if
               end if
               if (cell_c > 0.0_dp .and. s_cr > 0.0_dp) then
                  ey = cell_c/s_cr
                  tbx = tbx - dx*ey
                  tby = tby - dy*ey
                  tbz = tbz - dz*ey
               end if
            end do
            ! dcell_own/dR_c: (B) when c is the owner itself.
            if (ic == 1) then
               ox = tbx
               oy = tby
               oz = tbz
            end if
            vx = c*((ox - w_own*(tax + tbx))*inv_total)
            vy = c*((oy - w_own*(tay + tby))*inv_total)
            vz = c*((oz - w_own*(taz + tbz))*inv_total)
            vox = vox + vx
            voy = voy + vy
            voz = voz + vz
            !$acc atomic update
            grad(1, j) = grad(1, j) + vx
            !$acc atomic update
            grad(2, j) = grad(2, j) + vy
            !$acc atomic update
            grad(3, j) = grad(3, j) + vz
         end do
         ! The point moves with its owner: minus the sum, onto the owner.
         !$acc atomic update
         grad(1, own) = grad(1, own) - vox
         !$acc atomic update
         grad(2, own) = grad(2, own) - voy
         !$acc atomic update
         grad(3, own) = grad(3, own) - voz
      end if

      do ia = 1, n_active
         mark(g, active(g, ia)) = 0
      end do
   end subroutine partition_grad_point

   pure subroutine pair_sij(i, j, n_atoms, atom_coords, atom_r, unit_vec, shift, inv_distance, &
                            scheme, sval, d_own, d_oth)
      !! S(i,j): the fraction of the (i,j) pair's cutoff that atom i keeps,
      !! and its derivative with respect to R_i (`d_own`) and R_j (`d_oth`) --
      !! both one-sided, since S(i,j) depends on nothing else. Built from
      !! the SAME canonical pair (p, q) = (min(i,j), max(i,j)) regardless of
      !! which order (i,j) is asked in, so S(i,j) and S(j,i) = 1 - S(i,j)
      !! are bitwise consistent by construction.
      integer, intent(in) :: i, j, n_atoms
      real(dp), intent(in) :: atom_coords(3, n_atoms), atom_r(n_atoms), unit_vec(3, n_atoms)
      real(dp), intent(in) :: shift(n_atoms, n_atoms), inv_distance(n_atoms, n_atoms)
      integer, intent(in) :: scheme
      real(dp), intent(out) :: sval
      real(dp), intent(out) :: d_own(3), d_oth(3)
      real(dp) :: mu, nu, s, ds_dnu, dnu_dmu, invd
      real(dp) :: rij(3), ds_p(3), ds_q(3)
      integer :: p, q

      p = min(i, j)
      q = max(i, j)
      invd = inv_distance(p, q)
      rij = atom_coords(:, p) - atom_coords(:, q)

      mu = (atom_r(p) - atom_r(q))*invd
      nu = mu + shift(p, q)*(1.0_dp - mu*mu)

      if (scheme == PARTITION_BECKE) then
         s = becke_cutoff(nu)
         ds_dnu = becke_cutoff_derivative(nu)
      else
         s = stratmann_cutoff(nu)
         ds_dnu = stratmann_cutoff_derivative(nu)
      end if

      dnu_dmu = 1.0_dp - 2.0_dp*shift(p, q)*mu

      ! mu depends on R_p and R_q only: through the distances from the
      ! point, and through the internuclear distance in the denominator.
      ds_p = ds_dnu*dnu_dmu*(-unit_vec(:, p)*invd - (atom_r(p) - atom_r(q))*invd**3*rij)
      ds_q = ds_dnu*dnu_dmu*(unit_vec(:, q)*invd + (atom_r(p) - atom_r(q))*invd**3*rij)

      if (i < j) then
         sval = s
         d_own = ds_p
         d_oth = ds_q
      else
         sval = 1.0_dp - s
         d_own = -ds_q
         d_oth = -ds_p
      end if
   end subroutine pair_sij

   pure subroutine size_adjustment(atomic_numbers, adjust, shift)
      !! Pairwise shift applied to mu so unequal atoms get unequal cells
      integer, intent(in) :: atomic_numbers(:)
      integer, intent(in) :: adjust
      real(dp), intent(out) :: shift(:, :)

      real(dp), allocatable :: radius(:)
      integer :: n_atoms, i, j
      real(dp) :: chi, a

      n_atoms = size(atomic_numbers)
      allocate (radius(n_atoms))

      select case (adjust)
      case (ADJUST_TREUTLER)
         do i = 1, n_atoms
            radius(i) = sqrt(bragg_radius(atomic_numbers(i))) + TINY_RADIUS
         end do
      case (ADJUST_BECKE)
         do i = 1, n_atoms
            radius(i) = bragg_radius(atomic_numbers(i)) + TINY_RADIUS
         end do
      case default
         radius = 1.0_dp
      end select

      shift = 0.0_dp
      if (adjust == ADJUST_NONE) return

      do i = 1, n_atoms
         do j = 1, n_atoms
            if (i == j) cycle
            chi = radius(i)/radius(j)
            ! a = u/(u^2-1) with u = (chi-1)/(chi+1), which reduces to this.
            a = 0.25_dp*(1.0_dp/chi - chi)
            shift(i, j) = max(-MAX_ADJUST, min(MAX_ADJUST, a))
         end do
      end do
   end subroutine size_adjustment

   pure function becke_cutoff(nu) result(s)
      !! Becke's smooth cutoff: three iterations of p(x) = x(3 - x^2)/2
      !$acc routine seq
      real(dp), intent(in) :: nu
      real(dp) :: s
      real(dp) :: f

      f = nu
      f = 0.5_dp*f*(3.0_dp - f*f)
      f = 0.5_dp*f*(3.0_dp - f*f)
      f = 0.5_dp*f*(3.0_dp - f*f)
      s = 0.5_dp*(1.0_dp - f)
   end function becke_cutoff

   pure function stratmann_cutoff(mu) result(s)
      !! Stratmann's cutoff, exactly 1 or 0 outside |mu| >= a
      !$acc routine seq
      real(dp), intent(in) :: mu
      real(dp) :: s
      real(dp) :: z, z2

      if (mu <= -STRATMANN_A) then
         s = 1.0_dp
      else if (mu >= STRATMANN_A) then
         s = 0.0_dp
      else
         z = mu/STRATMANN_A
         z2 = z*z
         z = (1.0_dp/16.0_dp)*(z*(35.0_dp + z2*(-35.0_dp + z2*(21.0_dp - 5.0_dp*z2))))
         s = 0.5_dp*(1.0_dp - z)
      end if
   end function stratmann_cutoff

   pure function becke_cutoff_derivative(nu) result(ds)
      !! d/dnu of `becke_cutoff`
      !!
      !! Each iteration is p(x) = x(3 - x^2)/2 with p'(x) = 3(1 - x^2)/2, so the
      !! chain rule over the three of them is the product of p' at the value
      !! going into each.
      !$acc routine seq
      real(dp), intent(in) :: nu
      real(dp) :: ds
      real(dp) :: f0, f1, f2

      f0 = nu
      f1 = 0.5_dp*f0*(3.0_dp - f0*f0)
      f2 = 0.5_dp*f1*(3.0_dp - f1*f1)

      ds = -0.5_dp*(1.5_dp*(1.0_dp - f2*f2)) &
           *(1.5_dp*(1.0_dp - f1*f1)) &
           *(1.5_dp*(1.0_dp - f0*f0))
   end function becke_cutoff_derivative

   pure function stratmann_cutoff_derivative(mu) result(ds)
      !! d/dmu of `stratmann_cutoff`, zero where the cutoff is flat
      !$acc routine seq
      real(dp), intent(in) :: mu
      real(dp) :: ds
      real(dp) :: z, z2

      if (mu <= -STRATMANN_A .or. mu >= STRATMANN_A) then
         ds = 0.0_dp
      else
         z = mu/STRATMANN_A
         z2 = z*z
         ds = -0.5_dp*(1.0_dp/16.0_dp) &
              *(35.0_dp + z2*(-105.0_dp + z2*(105.0_dp - 35.0_dp*z2)))/STRATMANN_A
      end if
   end function stratmann_cutoff_derivative

end module trc_dft_partition
