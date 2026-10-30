!
! Hirshfeld atomic charges, by stockholder partition of the molecular density
!
! Ported from metalquicha's mqc_czt_hirshfeld, mqc_ -> trc_ renamed. The
! free-atom RECIPE and the renormalisation are mqc's, verbatim; the grid
! integral itself is written against terco's own grid/collocation
! machinery (trc_dft_grid, trc_collocation) rather than mqc's AO driver.
!
module trc_hirshfeld
   !! Each point of space is shared between the atoms in proportion to how
   !! much density each would put there on its own:
   !!
   !!     w_A(r) = rho_A^0(|r - R_A|) / sum_B rho_B^0(|r - R_B|)
   !!     q_A    = Z_A - integral w_A(r) rho(r) dr
   !!
   !! (Hirshfeld, Theor. Chim. Acta 44, 129 (1977)). The weights depend on
   !! the geometry and the elements only, never on the density being
   !! divided.
   !!
   !! **The free atoms are one fixed recipe, whatever the molecule ran.**
   !! Each element is the spherically averaged, spin-polarised atom of
   !! mqc's `mqc_czt_spherical_atom`, solved with wB97X in def2-TZVP on a
   !! grid fine enough to agree with PySCF's own 200x1202 to 2e-11 hartree
   !! (see mqc_czt_hirshfeld.f90's own header). terco does not solve this
   !! atom itself: `trc_hirshfeld_data` loads the already-converged radial
   !! density from test/basis_sets/hirshfeld_proatoms.json, generated
   !! offline by scripts/gen_hirshfeld_proatoms.py. A charge is therefore a
   !! property of the molecular density alone and does not move with the
   !! molecule's basis or method.
   !!
   !! The populations are renormalised to the density's exact electron
   !! count, Tr(D S): the grid misses a little of a diffuse density, and
   !! the charges then still sum to the molecular charge.
   !!
   !! terco has no effective-core-potential support at all (every basis it
   !! reads is all-electron), so mqc's own ECP refusal has nothing to
   !! check here. A ghost atom (zero nuclear charge) carries no free atom
   !! and takes no share, exactly as mqc's.
   use pic_types, only: dp
   use pic_logger, only: logger => global_logger
   use pic_strings, only: to_string
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_basis_types, only: trc_basis_t, trc_pairlist_t
   use trc_integrals, only: trc_1e
   use trc_dft_grid, only: dft_grid_t, dft_grid_build
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc, only: trc_xc_rho
   use trc_hirshfeld_data, only: trc_hirshfeld_table_t, trc_hirshfeld_ensure_loaded, HIRSHFELD_Z_MAX
   implicit none
   private

   public :: trc_hirshfeld_populations
   public :: trc_hirshfeld_free_atom_density
   public :: HIRSHFELD_GRID_LEVEL

   integer, parameter :: HIRSHFELD_GRID_LEVEL = 5
      !! Quadrature level for the molecular integral, when the caller names
      !! none -- mqc's own default (mqc_czt_hirshfeld.f90's
      !! HIRSHFELD_GRID_LEVEL), independent of whatever grid level the SCF
      !! used for its own XC integration.
   integer, parameter :: HIRSHFELD_BATCH = 256
      !! Points per spatial batch of the density evaluation.
   integer, parameter :: HIRSHFELD_CHUNK = 65536
      !! Points per stockholder-share pass: bounds share(natm, :) at natm*0.5 MB.

   real(dp), parameter :: GRID_DEFICIT_WARNING = 1.0e-4_dp
      !! How far the grid's raw electron count may miss the exact one
      !! before it is worth a warning rather than just the verbose log.

   type(trc_hirshfeld_table_t), save :: g_table
      !! The proatom tables, loaded once per process and kept for the life
      !! of it -- mqc's own `tables(Z_MAX)` module cache, the Fortran
      !! analogue.

contains

   subroutine trc_hirshfeld_populations(b, density, charges, error, spin_density, &
                                        spin_populations, grid_level, data_path)
      !! q_A = Z_A - N_A, with N_A the stockholder share of `density`
      !!
      !! `density` is the total, both spins, in terco's own CARTESIAN AO
      !! space (b%nao) -- the same convention `trc_xc_rks`/`trc_xc_uks`
      !! read (`scf_state%d_c`): a caller whose SCF ran in the spherical
      !! space (`b%spherical`) converts with `trc_sph`'s
      !! `sph_sph_to_cart` first (`trc_capi`'s `trc_hirshfeld_charges` does
      !! this). With `spin_density` (alpha minus beta, same space) the same
      !! weights give `spin_populations`.
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: density(b%nao, b%nao)
      real(dp), allocatable, intent(out) :: charges(:)   !! (natm), e
      type(error_t), intent(inout) :: error
      real(dp), intent(in), optional :: spin_density(b%nao, b%nao)
      real(dp), allocatable, intent(out), optional :: spin_populations(:)
      integer, intent(in), optional :: grid_level        !! Default `HIRSHFELD_GRID_LEVEL`
      character(len=*), intent(in), optional :: data_path
         !! Overrides the proatom table's default path

      type(dft_grid_t) :: grid
      type(trc_xc_grid_t) :: xg
      type(trc_basis_t) :: b_dev
      type(trc_pairlist_t) :: pl
      real(dp), allocatable :: smat(:, :), tmat(:, :), vmat(:, :)
      real(dp), allocatable :: rho(:), spin_rho(:), share(:, :)
      real(dp) :: wr
      integer :: p0, np_chunk, g
      real(dp), allocatable :: electrons(:), spins(:)
      real(dp), allocatable :: r_cut_list(:), log_rho_list(:, :), slope_list(:, :)
      logical, allocatable :: real_atom(:)
      real(dp) :: n_exact, spin_exact, deficit
      integer :: level, ia, z
      logical :: want_spin

      if (error%has_error()) return
      want_spin = present(spin_density) .and. present(spin_populations)

      call trc_hirshfeld_ensure_loaded(g_table, error, path=data_path)
      if (error%has_error()) then
         call error%set(error%get_code(), "hirshfeld charges: "//error%get_message())
         return
      end if

      allocate (real_atom(b%natm))
      do ia = 1, b%natm
         real_atom(ia) = b%at_z(ia) > 0.0_dp
         if (.not. real_atom(ia)) cycle
         z = b%at_elem(ia)
         if (z < 1 .or. z > HIRSHFELD_Z_MAX) then
            call error%set(ERROR_VALIDATION, "hirshfeld charges: no free atom for Z="// &
                           to_string(z)//"; the recipe covers H-Kr")
            return
         end if
         if (.not. g_table%valid(z)) then
            call error%set(ERROR_VALIDATION, "hirshfeld charges: the proatom table has no "// &
                           "entry for Z="//to_string(z))
            return
         end if
      end do

      allocate (r_cut_list(b%natm), log_rho_list(g_table%n_radial, b%natm), &
               slope_list(g_table%n_radial, b%natm))
      r_cut_list = 0.0_dp
      log_rho_list = 0.0_dp
      slope_list = 0.0_dp
      do ia = 1, b%natm
         if (.not. real_atom(ia)) cycle
         z = b%at_elem(ia)
         r_cut_list(ia) = g_table%r_cut(z)
         log_rho_list(:, ia) = g_table%log_rho(:, z)
         slope_list(:, ia) = g_table%slope(:, z)
      end do

      level = HIRSHFELD_GRID_LEVEL
      if (present(grid_level)) level = grid_level
      call dft_grid_build(b%at_r, b%at_elem, grid, error, level=level)
      if (error%has_error()) return

      ! rho through the XC integrator's batched path (screened per batch,
      ! one GEMM per batch shape), not a per-point nao^2 loop: on an
      ! 89-atom def2-TZVPD molecule that is ~1e11 operations instead of
      ! ~1e13, and no per-thread scratch the size of the basis. Every
      ! point is kept (w_tol 0) and functions are cut at 1e-14, far below
      ! what a charge can see. The grid comes back permuted; the weights
      ! and coordinates used below are xg's, in the same order as rho.
      call xg%build(grid%n_points, grid%coords, grid%weights, b, HIRSHFELD_BATCH, 1.0e-14_dp, 0.0_dp)
      call grid%destroy()
      call xg%to_device()
      if (b%on_device) then
         call density_on_grid(b, xg, density, rho)
         if (want_spin) call density_on_grid(b, xg, spin_density, spin_rho)
      else
         ! A caller outside an SCF context may hand over a host-only basis;
         ! the integrator reads it on the device.
         b_dev = b
         b_dev%on_device = .false.
         call b_dev%to_device()
         call density_on_grid(b_dev, xg, density, rho)
         if (want_spin) call density_on_grid(b_dev, xg, spin_density, spin_rho)
         call b_dev%release()
      end if

      ! The stockholder shares in chunks of points, so their (natm, points)
      ! storage stays bounded however large the grid.
      allocate (electrons(b%natm), spins(b%natm))
      electrons = 0.0_dp
      spins = 0.0_dp
      allocate (share(b%natm, HIRSHFELD_CHUNK))
      do p0 = 1, xg%npts, HIRSHFELD_CHUNK
         np_chunk = min(HIRSHFELD_CHUNK, xg%npts - p0 + 1)
         call hirshfeld_share_kernel(np_chunk, xg%r(:, p0:p0 + np_chunk - 1), b%natm, b%at_r, real_atom, &
                                     g_table%n_radial, g_table%log_r0, g_table%step, r_cut_list, &
                                     log_rho_list, slope_list, share)
         do g = 1, np_chunk
            wr = xg%w(p0 + g - 1)*rho(p0 + g - 1)
            do ia = 1, b%natm
               electrons(ia) = electrons(ia) + wr*share(ia, g)
            end do
         end do
         if (want_spin) then
            do g = 1, np_chunk
               wr = xg%w(p0 + g - 1)*spin_rho(p0 + g - 1)
               do ia = 1, b%natm
                  spins(ia) = spins(ia) + wr*share(ia, g)
               end do
            end do
         end if
      end do
      call xg%release()

      ! Sum to the density's own electron count, Tr(D S), not to what the
      ! grid integrates -- a diffuse density's tail is short by a little,
      ! and a difference of two states' charges (a condensed Fukui index,
      ! if one is ever built on this) inherits that and fails its sum rule
      ! otherwise.
      call pl%build(b, 0.0_dp)
      call pl%to_device()
      allocate (smat(b%nao, b%nao), tmat(b%nao, b%nao), vmat(b%nao, b%nao))
      call trc_1e(b, pl, smat, tmat, vmat)
      call pl%release()
      n_exact = sum(density*smat)
      if (want_spin) spin_exact = sum(spin_density*smat)
      call hirshfeld_renormalise(electrons, n_exact, deficit, want_spin, spins, spin_exact)
      call report_deficit(deficit)

      allocate (charges(b%natm))
      charges = 0.0_dp
      where (real_atom) charges = b%at_z - electrons
      if (want_spin) spin_populations = spins
   end subroutine trc_hirshfeld_populations

   subroutine hirshfeld_renormalise(electrons, n_exact, deficit, want_spin, spins, spin_exact)
      ! mqc_population_analysis's `renormalise_populations`, inlined: a
      ! ~10-line routine not worth a module of its own here, kept to
      ! exactly the same arithmetic.
      real(dp), intent(inout) :: electrons(:)
      real(dp), intent(in) :: n_exact
      real(dp), intent(out) :: deficit
      logical, intent(in) :: want_spin
      real(dp), intent(inout) :: spins(:)
      real(dp), intent(in) :: spin_exact

      real(dp) :: total, spin_deficit

      total = sum(electrons)
      deficit = n_exact - total
      if (total <= 0.0_dp) return
      if (want_spin) then
         spin_deficit = spin_exact - sum(spins)
         spins = spins + spin_deficit*electrons/total
      end if
      electrons = electrons*(n_exact/total)
   end subroutine hirshfeld_renormalise

   subroutine report_deficit(deficit)
      real(dp), intent(in) :: deficit

      call logger%verbose("  hirshfeld: the grid missed "//to_string(deficit)// &
                          " electrons, shared out by population")
      if (abs(deficit) > GRID_DEFICIT_WARNING) then
         call logger%warning("  hirshfeld: the integration grid missed "// &
                             to_string(deficit)//" electrons; the charges are "// &
                             "renormalised, but a grid this coarse for the density "// &
                             "is worth a finer one")
      end if
   end subroutine report_deficit

   subroutine trc_hirshfeld_free_atom_density(atomic_number, radii, rho, error, data_path)
      !! The free atom's total density at each radius, Bohr in, e/Bohr^3
      !! out -- for the tests, mqc's own `free_atom_density`.
      integer, intent(in) :: atomic_number
      real(dp), intent(in) :: radii(:)
      real(dp), allocatable, intent(out) :: rho(:)
      type(error_t), intent(inout) :: error
      character(len=*), intent(in), optional :: data_path

      integer :: k

      call trc_hirshfeld_ensure_loaded(g_table, error, path=data_path)
      if (error%has_error()) return
      if (atomic_number < 1 .or. atomic_number > HIRSHFELD_Z_MAX) then
         call error%set(ERROR_VALIDATION, "no free atom for Z="//to_string(atomic_number)// &
                        "; Hirshfeld covers H-Kr")
         return
      end if
      if (.not. g_table%valid(atomic_number)) then
         call error%set(ERROR_VALIDATION, "the proatom table has no entry for Z="// &
                        to_string(atomic_number))
         return
      end if

      allocate (rho(size(radii)))
      do k = 1, size(radii)
         rho(k) = hirshfeld_table_density(g_table%n_radial, g_table%log_r0, g_table%step, &
                                          g_table%r_cut(atomic_number), &
                                          g_table%log_rho(:, atomic_number), &
                                          g_table%slope(:, atomic_number), radii(k))
      end do
   end subroutine trc_hirshfeld_free_atom_density

   !
   ! KERNELS. Both are plain `do concurrent` over grid points, each thread
   ! calling a same-module `pure` `!$acc routine seq` item routine with
   ! explicit-shape dummies and its scratch declared inside -- the rule
   ! AGENTS.md states and `trc_xc.F90`/`trc_collocation.F90` follow
   ! throughout. Unscreened: every shell is evaluated at every point, which
   ! is the right trade for a one-shot post-SCF analysis on the molecule
   ! sizes Hirshfeld charges are asked for, not a per-iteration SCF kernel.
   !

   subroutine density_on_grid(b, xg, density, rho)
      !! rho at every point of `xg`, in its order; `b` and `xg` on the device.
      type(trc_basis_t), intent(in) :: b
      type(trc_xc_grid_t), intent(in) :: xg
      real(dp), intent(in) :: density(b%nao, b%nao)
      real(dp), allocatable, intent(out) :: rho(:)
      allocate (rho(xg%npts))
      call trc_xc_rho(b, xg, density, rho)
   end subroutine density_on_grid

   subroutine hirshfeld_share_kernel(n_points, coords, natm, at_r, real_atom, n_radial, log_r0, &
                                     step, r_cut, log_rho, slope, share)
      !! w_A at every grid point, for every atom -- zero for a ghost and
      !! zero at a point no free atom reaches.
      integer, intent(in) :: n_points, natm, n_radial
      real(dp), intent(in) :: coords(3, n_points), at_r(3, natm)
      logical, intent(in) :: real_atom(natm)
      real(dp), intent(in) :: log_r0, step
      real(dp), intent(in) :: r_cut(natm)
      real(dp), intent(in) :: log_rho(n_radial, natm), slope(n_radial, natm)
      real(dp), intent(out) :: share(natm, n_points)
      integer :: g

      !$acc enter data copyin(coords, at_r, real_atom, r_cut, log_rho, slope) create(share)
      do concurrent(g=1:n_points)
         call hirshfeld_share_point(g, n_points, coords, natm, at_r, real_atom, n_radial, log_r0, &
                                    step, r_cut, log_rho, slope, share(:, g))
      end do
      !$acc update self(share)
      !$acc exit data delete(coords, at_r, real_atom, r_cut, log_rho, slope, share)
   end subroutine hirshfeld_share_kernel

   pure subroutine hirshfeld_share_point(g, n_points, coords, natm, at_r, real_atom, n_radial, &
                                         log_r0, step, r_cut, log_rho, slope, share_col)
      !$acc routine seq
      integer, intent(in) :: g, n_points, natm, n_radial
      real(dp), intent(in) :: coords(3, n_points), at_r(3, natm)
      logical, intent(in) :: real_atom(natm)
      real(dp), intent(in) :: log_r0, step
      real(dp), intent(in) :: r_cut(natm)
      real(dp), intent(in) :: log_rho(n_radial, natm), slope(n_radial, natm)
      real(dp), intent(out) :: share_col(natm)

      real(dp) :: total
      real(dp) :: d(3)
      real(dp) :: r
      integer :: ia

      ! The free-atom densities go straight into share_col, then are
      ! divided: no per-thread scratch of length natm on the device.
      total = 0.0_dp
      do ia = 1, natm
         share_col(ia) = 0.0_dp
         if (.not. real_atom(ia)) cycle
         d(1) = coords(1, g) - at_r(1, ia)
         d(2) = coords(2, g) - at_r(2, ia)
         d(3) = coords(3, g) - at_r(3, ia)
         r = sqrt(d(1)*d(1) + d(2)*d(2) + d(3)*d(3))
         share_col(ia) = hirshfeld_table_density(n_radial, log_r0, step, r_cut(ia), &
                                                 log_rho(:, ia), slope(:, ia), r)
         total = total + share_col(ia)
      end do
      if (total > 0.0_dp) then
         do ia = 1, natm
            share_col(ia) = share_col(ia)/total
         end do
      end if
   end subroutine hirshfeld_share_point

   pure function hirshfeld_table_density(n_radial, log_r0, step, r_cut, log_rho_col, slope_col, r) &
      result(rho)
      !! Cubic Hermite in (ln r, ln rho) between the two bracketing radii --
      !! mqc_czt_hirshfeld.f90's `table_density`, verbatim arithmetic.
      !$acc routine seq
      integer, intent(in) :: n_radial
      real(dp), intent(in) :: log_r0, step, r_cut, r
      real(dp), intent(in) :: log_rho_col(n_radial), slope_col(n_radial)
      real(dp) :: rho

      real(dp) :: x, t, h00, h10, h01, h11
      integer :: k

      rho = 0.0_dp
      if (r >= r_cut) return
      x = (log(max(r, exp(log_r0))) - log_r0)/step
      k = min(int(x) + 1, n_radial - 1)
      t = x - real(k - 1, dp)
      h00 = (1.0_dp + 2.0_dp*t)*(1.0_dp - t)**2
      h10 = t*(1.0_dp - t)**2
      h01 = t*t*(3.0_dp - 2.0_dp*t)
      h11 = t*t*(t - 1.0_dp)
      rho = exp(h00*log_rho_col(k) + h10*step*slope_col(k) + &
                h01*log_rho_col(k + 1) + h11*step*slope_col(k + 1))
   end function hirshfeld_table_density

end module trc_hirshfeld
