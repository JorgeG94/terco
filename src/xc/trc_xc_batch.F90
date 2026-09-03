!
! The grid, batched and screened for the XC integrator.
!
! A molecular grid is tens of thousands of points, and a basis function has
! support on a small fraction of them. Evaluating every function at every
! point is the arithmetic the device does faster without getting anywhere:
! it scales as n_points x n_ao when the useful work scales as n_points x
! (functions nearby). This container is the difference between the two.
!
! The design is GauXC's (BSD; Williams-Young et al.), taken as a design and
! not as code:
!
!   1. cut the points into compact spatial BATCHES, by splitting the
!      bounding box along its longest axis at its midpoint until a box holds
!      at most `max_pts` points;
!   2. give every shell a CUTOFF RADIUS beyond which all its primitives are
!      below `tol`;
!   3. for each batch keep only the shells whose radius reaches its box.
!
! Everything downstream then works on a batch's LOCAL basis: `b_sh` lists
! its shells, `b_ao` the global AO index of each local function, and the
! kernels never see the rest.
!
! The points are stored PERMUTED, batch after batch, so that a batch is a
! contiguous range and a thread's neighbours in the launch are its
! neighbours in space. `batch_of` maps a point back to its batch.
!
module trc_xc_batch
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t
   implicit none
   private

   public :: trc_xc_grid_t

   type :: trc_xc_grid_t
      integer :: npts = 0, nbatch = 0, max_nloc = 0
      real(dp), allocatable :: r(:, :)      !! (3, npts), permuted
      real(dp), allocatable :: w(:)         !! (npts), permuted
      integer, allocatable :: batch_of(:)   !! (npts)
      integer, allocatable :: b_off(:)      !! (nbatch+1): batch b is points b_off(b) .. b_off(b+1)-1
      integer, allocatable :: b_shoff(:)    !! (nbatch+1): its shells are b_sh(b_shoff(b) .. b_shoff(b+1)-1)
      integer, allocatable :: b_sh(:)
      integer, allocatable :: b_aooff(:)    !! (nbatch+1): its AOs are b_ao(b_aooff(b) .. b_aooff(b+1)-1)
      integer, allocatable :: b_ao(:)
      logical :: on_device = .false.
   contains
      procedure :: build => xcgrid_build
      procedure :: to_device => xcgrid_to_device
      procedure :: release => xcgrid_release
   end type trc_xc_grid_t

contains

   pure integer function ncart(l)
      integer, intent(in) :: l
      ncart = (l + 1)*(l + 2)/2
   end function ncart

   !
   ! Radius beyond which |c| r^l exp(-a r^2) < tol, for the primitive that
   ! reaches furthest. Found by bisection outward from the maximum of
   ! r^l exp(-a r^2), which sits at r^2 = l/(2a); a Gaussian is monotone
   ! past it, so the bracket is sound.
   !
   pure function shell_cutoff(l, np, e, c, tol) result(rcut)
      integer, intent(in) :: l, np
      real(dp), intent(in) :: e(np), c(np), tol
      real(dp) :: rcut
      real(dp) :: lo, hi, mid, f, ac
      integer :: p, it

      rcut = 0.0_dp
      do p = 1, np
         ac = abs(c(p))
         if (ac <= 0.0_dp) cycle
         lo = sqrt(real(l, dp)/(2.0_dp*e(p)))
         if (envelope(lo) < tol) cycle
         hi = lo + 1.0_dp
         do while (envelope(hi) >= tol)
            hi = 2.0_dp*hi + 1.0_dp
         end do
         do it = 1, 60
            mid = 0.5_dp*(lo + hi)
            f = envelope(mid)
            if (f >= tol) then
               lo = mid
            else
               hi = mid
            end if
            if (hi - lo < 1.0e-6_dp) exit
         end do
         rcut = max(rcut, hi)
      end do
   contains
      pure real(dp) function envelope(r)
         real(dp), intent(in) :: r
         envelope = ac*r**l*exp(-e(p)*r*r)
      end function envelope
   end function shell_cutoff

   !
   ! Build from a grid and a basis.
   !
   ! `tol` is the value below which a basis function is treated as zero on a
   ! batch. 1e-10 is what GauXC ships; the error it introduces in an energy
   ! is far below the grid's own, and check_xc_energy measures it against
   ! the unscreened evaluation rather than assuming.
   !
   subroutine xcgrid_build(this, npts, coords, weights, b, max_pts, tol)
      class(trc_xc_grid_t), intent(inout) :: this
      integer, intent(in) :: npts
      real(dp), intent(in) :: coords(3, npts), weights(npts)
      type(trc_basis_t), intent(in) :: b
      integer, intent(in), optional :: max_pts
      real(dp), intent(in), optional :: tol

      integer :: maxp, i, j, k, lo, hi, axis, top, nb, ib, ish, nsh_b, nao_b
      integer :: sh_in_batch
      integer, allocatable :: idx(:), st_lo(:), st_hi(:), b_lo(:), b_hi(:)
      real(dp), allocatable :: rcut(:), bmin(:, :), bmax(:, :)
      real(dp) :: lo3(3), hi3(3), mid, d2, dd
      real(dp) :: eps_tol

      call this%release()
      maxp = 512
      if (present(max_pts)) maxp = max(1, max_pts)
      eps_tol = 1.0e-10_dp
      if (present(tol)) eps_tol = tol
      this%npts = npts
      if (npts == 0) return

      ! --- 1. spatial batches, by recursive midpoint bisection ---------------
      allocate (idx(npts), st_lo(64), st_hi(64))
      allocate (b_lo(npts), b_hi(npts))   ! at most npts batches
      do i = 1, npts
         idx(i) = i
      end do
      nb = 0
      top = 1
      st_lo(1) = 1; st_hi(1) = npts
      do while (top > 0)
         lo = st_lo(top); hi = st_hi(top); top = top - 1
         if (hi - lo + 1 <= maxp) then
            nb = nb + 1
            b_lo(nb) = lo; b_hi(nb) = hi
            cycle
         end if
         lo3 = huge(1.0_dp); hi3 = -huge(1.0_dp)
         do k = lo, hi
            lo3 = min(lo3, coords(:, idx(k)))
            hi3 = max(hi3, coords(:, idx(k)))
         end do
         axis = maxloc(hi3 - lo3, dim=1)
         mid = 0.5_dp*(lo3(axis) + hi3(axis))
         ! Partition idx(lo:hi): below the midpoint first.
         i = lo; j = hi
         do while (i <= j)
            if (coords(axis, idx(i)) < mid) then
               i = i + 1
            else
               k = idx(i); idx(i) = idx(j); idx(j) = k
               j = j - 1
            end if
         end do
         ! i is the first index at or above the midpoint. A box of coincident
         ! points cannot be split by value; split it by count instead.
         if (i <= lo .or. i > hi) i = (lo + hi)/2 + 1
         if (top + 2 > size(st_lo)) call grow_stack()
         top = top + 1; st_lo(top) = lo; st_hi(top) = i - 1
         top = top + 1; st_lo(top) = i; st_hi(top) = hi
      end do
      this%nbatch = nb

      allocate (this%r(3, npts), this%w(npts), this%batch_of(npts), this%b_off(nb + 1))
      k = 0
      do ib = 1, nb
         this%b_off(ib) = k + 1
         do i = b_lo(ib), b_hi(ib)
            k = k + 1
            this%r(:, k) = coords(:, idx(i))
            this%w(k) = weights(idx(i))
            this%batch_of(k) = ib
         end do
      end do
      this%b_off(nb + 1) = npts + 1

      ! --- 2. shell cutoff radii --------------------------------------------
      allocate (rcut(b%nshell))
      do ish = 1, b%nshell
         rcut(ish) = shell_cutoff(b%sh_l(ish), b%sh_np(ish), b%sh_e(1:b%sh_np(ish), ish), &
                                  b%sh_c(1:b%sh_np(ish), ish), eps_tol)
      end do

      ! --- 3. the local basis of every batch ---------------------------------
      allocate (bmin(3, nb), bmax(3, nb))
      do ib = 1, nb
         bmin(:, ib) = huge(1.0_dp); bmax(:, ib) = -huge(1.0_dp)
         do k = this%b_off(ib), this%b_off(ib + 1) - 1
            bmin(:, ib) = min(bmin(:, ib), this%r(:, k))
            bmax(:, ib) = max(bmax(:, ib), this%r(:, k))
         end do
      end do
      allocate (this%b_shoff(nb + 1), this%b_aooff(nb + 1))
      ! Count, then fill.
      this%b_shoff(1) = 1; this%b_aooff(1) = 1
      this%max_nloc = 0
      do ib = 1, nb
         nsh_b = 0; nao_b = 0
         do ish = 1, b%nshell
            if (reaches(ish, ib)) then
               nsh_b = nsh_b + 1
               nao_b = nao_b + ncart(b%sh_l(ish))
            end if
         end do
         this%b_shoff(ib + 1) = this%b_shoff(ib) + nsh_b
         this%b_aooff(ib + 1) = this%b_aooff(ib) + nao_b
         this%max_nloc = max(this%max_nloc, nao_b)
      end do
      allocate (this%b_sh(this%b_shoff(nb + 1) - 1), this%b_ao(this%b_aooff(nb + 1) - 1))
      do ib = 1, nb
         sh_in_batch = this%b_shoff(ib)
         k = this%b_aooff(ib)
         do ish = 1, b%nshell
            if (reaches(ish, ib)) then
               this%b_sh(sh_in_batch) = ish
               sh_in_batch = sh_in_batch + 1
               do i = 0, ncart(b%sh_l(ish)) - 1
                  this%b_ao(k) = b%sh_ao(ish) + i
                  k = k + 1
               end do
            end if
         end do
      end do

   contains

      logical function reaches(ish, ib)
         integer, intent(in) :: ish, ib
         integer :: dim
         d2 = 0.0_dp
         do dim = 1, 3
            dd = max(0.0_dp, bmin(dim, ib) - b%sh_r(dim, ish), b%sh_r(dim, ish) - bmax(dim, ib))
            d2 = d2 + dd*dd
         end do
         reaches = d2 < rcut(ish)*rcut(ish)
      end function reaches

      subroutine grow_stack()
         integer, allocatable :: t(:)
         allocate (t(2*size(st_lo)))
         t(1:size(st_lo)) = st_lo
         call move_alloc(t, st_lo)
         allocate (t(2*size(st_hi)))
         t(1:size(st_hi)) = st_hi
         call move_alloc(t, st_hi)
      end subroutine grow_stack

   end subroutine xcgrid_build

   subroutine xcgrid_to_device(this)
      class(trc_xc_grid_t), intent(inout) :: this
      if (this%on_device .or. this%npts == 0) return
      !$acc enter data copyin(this)
      !$acc enter data copyin(this%r, this%w, this%batch_of, this%b_off, &
      !$acc                   this%b_shoff, this%b_sh, this%b_aooff, this%b_ao)
      this%on_device = .true.
   end subroutine xcgrid_to_device

   subroutine xcgrid_release(this)
      class(trc_xc_grid_t), intent(inout) :: this
      if (this%on_device) then
         !$acc exit data delete(this%r, this%w, this%batch_of, this%b_off, &
         !$acc                  this%b_shoff, this%b_sh, this%b_aooff, this%b_ao)
         !$acc exit data delete(this)
         this%on_device = .false.
      end if
      if (allocated(this%r)) deallocate (this%r)
      if (allocated(this%w)) deallocate (this%w)
      if (allocated(this%batch_of)) deallocate (this%batch_of)
      if (allocated(this%b_off)) deallocate (this%b_off)
      if (allocated(this%b_shoff)) deallocate (this%b_shoff)
      if (allocated(this%b_sh)) deallocate (this%b_sh)
      if (allocated(this%b_aooff)) deallocate (this%b_aooff)
      if (allocated(this%b_ao)) deallocate (this%b_ao)
      this%npts = 0; this%nbatch = 0; this%max_nloc = 0
   end subroutine xcgrid_release

end module trc_xc_batch
