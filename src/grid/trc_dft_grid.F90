module trc_dft_grid
   !! Molecular integration grid for DFT
   !!
   !! Assembles the three pieces into one quadrature. Each atom gets a spherical
   !! product grid -- a Treutler-Ahlrichs radial mesh times a Lebedev sphere --
   !! and the Becke partition then decides how much of each point belongs to the
   !! atom that produced it, so the atomic grids together cover space once.
   !!
   !! The weight of a point is
   !!
   !!    w = 4*pi * r^2 * dr * w_lebedev * w_becke
   !!
   !! where w_lebedev sums to 1 over the sphere. The 4*pi therefore appears
   !! exactly once, here; `trc_dft_radial` deliberately leaves it out of its
   !! weights so it cannot be applied twice.
   !!
   !! Grid size follows a `level` from 0 to 9, with radial and angular counts
   !! chosen per element period. These are the same tables PySCF uses, which is
   !! what makes a level-for-level comparison meaningful. Level 3 is the usual
   !! default and is the one a production calculation should start from.
   use pic_types, only: dp, int64
   use pic_timer, only: timer_type
   use trc_constants, only: PI
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_lebedev, only: lebedev_grid
   use trc_dft_radial, only: dft_radial_treutler_ahlrichs
   use trc_dft_partition, only: dft_partition_becke_weights, dft_partition_options_t
   use trc_dft_prune, only: dft_prune_angular_orders, PRUNE_NONE, PRUNE_NWCHEM
   implicit none
   private

   public :: dft_grid_t, dft_grid_options_t
   public :: dft_grid_build, dft_grid_build_block
   public :: dft_grid_level_radial, dft_grid_level_angular
   public :: MIN_GRID_LEVEL, MAX_GRID_LEVEL, DEFAULT_GRID_LEVEL
   public :: DEFAULT_PRUNE

   integer, parameter :: N_DIM = 3
   integer, parameter :: N_PERIODS = 7
   integer, parameter :: MIN_GRID_LEVEL = 0
   integer, parameter :: MAX_GRID_LEVEL = 9
   integer, parameter :: DEFAULT_GRID_LEVEL = 3

   integer, parameter :: DEFAULT_PRUNE = PRUNE_NWCHEM
   !! Pruning is on by default, as it is in every production DFT code. An
   !! unpruned level-3 water grid is 52850 points against 20240 pruned, for the
   !! same accuracy to well past where the functional itself is reliable.

   integer, parameter :: PERIOD_LAST_Z(N_PERIODS) = [2, 10, 18, 36, 54, 86, 118]
   !! Highest Z in each period, used to place an element in the level tables

   integer, parameter :: RAD_GRIDS(MIN_GRID_LEVEL:MAX_GRID_LEVEL, N_PERIODS) = &
                         transpose(reshape([ &
                                           10, 15, 20, 30, 35, 40, 50, &
                                           30, 40, 50, 60, 65, 70, 75, &
                                           40, 60, 65, 75, 80, 85, 90, &
                                           50, 75, 80, 90, 95, 100, 105, &
                                           60, 90, 95, 105, 110, 115, 120, &
                                           70, 105, 110, 120, 125, 130, 135, &
                                           80, 120, 125, 135, 140, 145, 150, &
                                           90, 135, 140, 150, 155, 160, 165, &
                                           100, 150, 155, 165, 170, 175, 180, &
                                           200, 200, 200, 200, 200, 200, 200], [N_PERIODS, MAX_GRID_LEVEL + 1]))
   !! Radial shells per (level, period)

   integer, parameter :: ANG_POINTS(MIN_GRID_LEVEL:MAX_GRID_LEVEL, N_PERIODS) = &
                         transpose(reshape([ &
                                           50, 86, 110, 110, 110, 110, 110, &
                                           110, 194, 194, 194, 194, 194, 194, &
                                           194, 302, 302, 302, 302, 302, 302, &
                                           302, 302, 434, 434, 434, 434, 434, &
                                           434, 590, 590, 590, 590, 590, 590, &
                                           590, 770, 770, 770, 770, 770, 770, &
                                           770, 974, 974, 974, 974, 974, 974, &
                                           974, 1202, 1202, 1202, 1202, 1202, 1202, &
                                           1202, 1202, 1202, 1202, 1202, 1202, 1202, &
                                           1454, 1454, 1454, 1454, 1454, 1454, 1454], &
                                           [N_PERIODS, MAX_GRID_LEVEL + 1]))
   !! Lebedev points per (level, period)
   !!
   !! PySCF's table stores Lebedev degrees and maps them to point counts; these
   !! are the point counts directly, since that is what `lebedev_grid` takes.

   type :: dft_grid_options_t
      !! How the grid is built, beyond its `level`. The defaults are the
      !! production grid.
      integer :: prune = DEFAULT_PRUNE
      !! PRUNE_NONE or PRUNE_NWCHEM.
      integer :: n_radial = 0, n_angular = 0
      !! Overrides for every atom, which is what a convergence study wants;
      !! zero leaves the level's tables in charge. Supplying one without the
      !! other is refused rather than silently half-applied.
      type(dft_partition_options_t) :: partition
      !! The partition's cutoff and size adjustment (and screening, scratch
      !! budget); its point block is set by dft_grid_build_block.
   end type dft_grid_options_t

   type :: dft_grid_t
      !! A molecular quadrature: points, weights, and which atom made each point
      integer :: n_points = 0
      real(dp), allocatable :: coords(:, :)  !! (3, n_points), Bohr
      real(dp), allocatable :: weights(:)    !! (n_points), full volume weight
      real(dp), allocatable :: quad_weights(:)
         !! (n_points), the radial times angular weight *before* the Becke
         !! partition multiplies it in. Kept because a gradient needs the
         !! partition weight and the quadrature weight separately -- the
         !! partition depends on where the nuclei are and the quadrature does
         !! not, so only one of the two carries a derivative. Recovering it by
         !! dividing `weights` by the partition would divide by zero exactly
         !! where the partition underflows.
      integer, allocatable :: atom(:)        !! (n_points), owning atom
      ! What the partition was built with. Recorded because a gradient has to
      ! differentiate the same partition the energy integrated over, and a
      ! caller reconstructing these from defaults would silently differentiate
      ! a different one the moment the defaults changed.
      type(dft_partition_options_t) :: partition
      !! Cutoff, size adjustment and screening; the point block is reset to
      !! every point, since a gradient walks all of them.
      integer, allocatable :: numbers(:)     !! (n_atoms), Z, for the partition radii
      real(dp) :: t_points = 0.0_dp, t_partition = 0.0_dp
      !! Wall time of the point-generation loop (radial x Lebedev x pruning)
      !! and of `apply_partition` (the Becke/Stratmann weights), for the
      !! caller's own "grid (s)" sub-table -- trc_scf_context.F90's `setup`.
   contains
      procedure :: destroy => dft_grid_destroy
   end type dft_grid_t

contains

   pure function element_period(atomic_number) result(period)
      !! Period of an element, 1 to 7, for indexing the level tables
      integer, intent(in) :: atomic_number
      integer :: period
      integer :: i

      period = N_PERIODS
      do i = 1, N_PERIODS
         if (atomic_number <= PERIOD_LAST_Z(i)) then
            period = i
            exit
         end if
      end do
   end function element_period

   pure function dft_grid_level_radial(atomic_number, level) result(n_radial)
      !! Radial shells this element gets at this level
      integer, intent(in) :: atomic_number, level
      integer :: n_radial
      n_radial = RAD_GRIDS(clamp_level(level), element_period(atomic_number))
   end function dft_grid_level_radial

   pure function dft_grid_level_angular(atomic_number, level) result(n_angular)
      !! Lebedev points this element gets at this level
      integer, intent(in) :: atomic_number, level
      integer :: n_angular
      n_angular = ANG_POINTS(clamp_level(level), element_period(atomic_number))
   end function dft_grid_level_angular

   pure function clamp_level(level) result(clamped)
      integer, intent(in) :: level
      integer :: clamped
      clamped = max(MIN_GRID_LEVEL, min(MAX_GRID_LEVEL, level))
   end function clamp_level

   subroutine dft_grid_build(atom_coords, atomic_numbers, grid, error, level, options)
      !! Build the molecular grid: see `build_grid` for the arguments.
      !!
      !! A separate entry from `dft_grid_build_block`, with this interface
      !! unchanged, because nvfortran 26.5 crashes compiling a caller that
      !! passes an expression actual with keyword optionals to a version
      !! carrying the two rank arguments (check_xc_energy did).
      real(dp), intent(in) :: atom_coords(:, :)
      integer, intent(in) :: atomic_numbers(:)
      type(dft_grid_t), intent(out) :: grid
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: level
      type(dft_grid_options_t), intent(in), optional :: options
      call build_grid(atom_coords, atomic_numbers, grid, error, level, options)
   end subroutine dft_grid_build

   subroutine dft_grid_build_block(atom_coords, atomic_numbers, grid, error, my_rank, n_ranks, level, options)
      !! Build the molecular grid with this rank partitioning only its block
      !! of the points: the partition is the expensive part of the build,
      !! and the same on every rank. The returned `weights` are zero outside
      !! the block, so the caller sums them across ranks; the quadrature
      !! weights are complete regardless.
      real(dp), intent(in) :: atom_coords(:, :)
      integer, intent(in) :: atomic_numbers(:)
      type(dft_grid_t), intent(out) :: grid
      type(error_t), intent(inout) :: error
      integer, intent(in) :: my_rank, n_ranks
      integer, intent(in), optional :: level
      type(dft_grid_options_t), intent(in), optional :: options
      call build_grid(atom_coords, atomic_numbers, grid, error, level, options, my_rank, n_ranks)
   end subroutine dft_grid_build_block

   subroutine build_grid(atom_coords, atomic_numbers, grid, error, level, options, my_rank, n_ranks)
      !! Build the molecular grid
      !!
      !! `level` picks per-element sizes from the standard tables;
      !! `options` overrides the pruning, the per-atom sizes and the
      !! partition.
      real(dp), intent(in) :: atom_coords(:, :)   !! (3, n_atoms), Bohr
      integer, intent(in) :: atomic_numbers(:)    !! Z per atom
      type(dft_grid_t), intent(out) :: grid
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: level      !! 0 to 9, default 3
      type(dft_grid_options_t), intent(in), optional :: options
      !! Default: the production grid (see dft_grid_options_t).
      integer, intent(in), optional :: my_rank, n_ranks
      !! With both given, only this rank's block of the points is partitioned.

      real(dp), allocatable :: r(:), dr(:), sphere(:, :), w_ang(:)
      integer, allocatable :: shell_order(:)
      type(dft_grid_options_t) :: opts
      integer :: n_atoms, used_level
      integer :: ia, i, j, k, nr, na, total, offset, cached_order
      type(timer_type) :: tm_gen

      n_atoms = size(atom_coords, 2)
      if (size(atomic_numbers) /= n_atoms) then
         call error%set(ERROR_VALIDATION, "grid: atomic_numbers does not match atom_coords")
         return
      end if
      if (n_atoms < 1) then
         call error%set(ERROR_VALIDATION, "grid: no atoms")
         return
      end if

      used_level = DEFAULT_GRID_LEVEL
      if (present(level)) used_level = level
      if (present(options)) opts = options

      if ((opts%n_radial > 0) .neqv. (opts%n_angular > 0)) then
         call error%set(ERROR_VALIDATION, &
                        "grid: n_radial and n_angular must be given together")
         return
      end if

      call tm_gen%start()
      ! Size the whole grid first so it can be allocated once. With pruning the
      ! count is no longer nr*na: each shell carries its own order, and those
      ! depend on the radii, so the radial mesh has to be built here as well.
      total = 0
      do ia = 1, n_atoms
         call atom_sizes(atomic_numbers(ia), used_level, opts, nr, na)
         call dft_radial_treutler_ahlrichs(nr, atomic_numbers(ia), r, dr, error)
         if (error%has_error()) return
         allocate (shell_order(nr))
         call dft_prune_angular_orders(opts%prune, atomic_numbers(ia), r, na, shell_order, error)
         if (error%has_error()) return
         total = total + sum(shell_order)
         deallocate (r, dr, shell_order)
      end do

      grid%n_points = total
      allocate (grid%coords(N_DIM, total))
      allocate (grid%weights(total))
      allocate (grid%atom(total))

      offset = 0
      do ia = 1, n_atoms
         call atom_sizes(atomic_numbers(ia), used_level, opts, nr, na)

         call dft_radial_treutler_ahlrichs(nr, atomic_numbers(ia), r, dr, error)
         if (error%has_error()) return

         allocate (shell_order(nr))
         call dft_prune_angular_orders(opts%prune, atomic_numbers(ia), r, na, shell_order, error)
         if (error%has_error()) return

         ! Shells sharing an order are not contiguous -- the outer zone repeats
         ! the order of the third -- so hold the last sphere built and reuse it
         ! when the order does not change, rather than regenerating a 302-point
         ! grid for every shell that wants one.
         cached_order = 0
         do i = 1, nr
            if (shell_order(i) /= cached_order) then
               if (allocated(sphere)) deallocate (sphere, w_ang)
               call lebedev_grid(shell_order(i), sphere, w_ang, error)
               if (error%has_error()) return
               cached_order = shell_order(i)
            end if

            do j = 1, shell_order(i)
               k = offset + j
               grid%coords(:, k) = atom_coords(:, ia) + r(i)*sphere(:, j)
               ! 4*pi appears once, here: w_ang sums to 1 and dr excludes it.
               grid%weights(k) = 4.0_dp*PI*r(i)*r(i)*dr(i)*w_ang(j)
               grid%atom(k) = ia
            end do
            offset = offset + shell_order(i)
         end do

         deallocate (r, dr, shell_order)
         if (allocated(sphere)) deallocate (sphere, w_ang)
      end do

      grid%partition = opts%partition
      grid%partition%first = 1
      grid%partition%last = huge(1)
      allocate (grid%numbers(n_atoms), source=atomic_numbers)
      grid%t_points = tm_gen%get_elapsed_time()

      call tm_gen%start()
      call apply_partition(grid, atom_coords, atomic_numbers, error, my_rank, n_ranks)
      grid%t_partition = tm_gen%get_elapsed_time()
   end subroutine build_grid

   pure subroutine atom_sizes(atomic_number, level, options, nr, na)
      !! Radial and angular counts for one atom, honouring any override
      integer, intent(in) :: atomic_number, level
      type(dft_grid_options_t), intent(in) :: options
      integer, intent(out) :: nr, na

      if (options%n_radial > 0) then
         nr = options%n_radial
         na = options%n_angular
      else
         nr = dft_grid_level_radial(atomic_number, level)
         na = dft_grid_level_angular(atomic_number, level)
      end if
   end subroutine atom_sizes

   subroutine apply_partition(grid, atom_coords, atomic_numbers, error, my_rank, n_ranks)
      !! Multiply the product weights by the Becke cell weights, built as
      !! grid%partition says
      type(dft_grid_t), intent(inout) :: grid
      real(dp), intent(in) :: atom_coords(:, :)
      integer, intent(in) :: atomic_numbers(:)
      type(error_t), intent(inout) :: error
      integer, intent(in), optional :: my_rank, n_ranks

      real(dp), allocatable :: cell(:)
      type(dft_partition_options_t) :: block_partition

      allocate (cell(grid%n_points))
      block_partition = grid%partition
      if (present(my_rank) .and. present(n_ranks)) then
         ! Contiguous blocks: the chunks inside are contiguous too.
         block_partition%first = int((int(grid%n_points, int64)*my_rank)/n_ranks) + 1
         block_partition%last = int((int(grid%n_points, int64)*(my_rank + 1))/n_ranks)
      end if
      call dft_partition_becke_weights(grid%coords, atom_coords, atomic_numbers, &
                                       grid%atom, block_partition, cell, error)
      if (error%has_error()) return

      ! Before, not after: this is the quadrature weight on its own.
      if (allocated(grid%quad_weights)) deallocate (grid%quad_weights)
      allocate (grid%quad_weights(grid%n_points), source=grid%weights)

      grid%weights = grid%weights*cell
   end subroutine apply_partition

   pure subroutine dft_grid_destroy(this)
      !! Release the grid
      class(dft_grid_t), intent(inout) :: this

      this%n_points = 0
      if (allocated(this%coords)) deallocate (this%coords)
      if (allocated(this%weights)) deallocate (this%weights)
      if (allocated(this%quad_weights)) deallocate (this%quad_weights)
      if (allocated(this%numbers)) deallocate (this%numbers)
      if (allocated(this%atom)) deallocate (this%atom)
   end subroutine dft_grid_destroy

end module trc_dft_grid
