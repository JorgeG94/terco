module check_xc_grad_response_suite
   !!
   !! trc_xc_grad_response_many (src/xc/trc_xc_grad_response.F90): for each
   !! column I, the nuclear derivative at fixed AO matrices of
   !!
   !!     G_I = Tr(P_I V_xc[D]) + 1/2 E_xc''[T_I, T_I],
   !!
   !! LDA and GGA, restricted, all columns in one sweep of the grid.
   !!
   !! Convention checked here, what the TDDFT gradient will rely on: the first
   !! term is Tr(P V_xc[D]) as trc_xc_rks forms V_xc, the second is
   !! 1/2 Tr(T f_xc[T]) with f_xc[T] as trc_xc_fxc_rks forms it (the
   !! second derivative with respect to the TOTAL density), so a caller who
   !! solves with trc_xc_fxc_rks_apply_many on a vector T gets the matching
   !! gradient from the same T, and the spin factor is the caller's (the term
   !! is quadratic in T). P and T enter through their symmetric parts.
   !!
   !! The gradient includes the grid response (the points move with their
   !! atoms and the Becke weights follow), as trc_xc_grad and
   !! xc_grad_potential_add do, so the finite differences REBUILD the grid at
   !! every displaced geometry, with the AO matrices D, P_I, T_I held fixed.
   !!
   !! Checks:
   !!  1. central finite differences in the nuclear coordinates (Richardson
   !!     over steps 2h and h) of Tr(P V_xc[D]) from trc_xc_rks and of
   !!     1/2 Tr(T f_xc[T]) from trc_xc_fxc_rks, rebuilt at each displaced
   !!     geometry: each term alone and both together, water and H2O2, in
   !!     6-31G*; plus translational invariance, and the negative controls
   !!     (without the k_xc part, without the grid response);
   !!  2. the identity in a density step, 1/2 d/dR E''[T,T] =
   !!     1/2 d^2/dt^2 of trc_xc_grad's gradient at D + t T (and the first
   !!     derivative for P), which needs no nuclear finite difference;
   !!  3. a batch of columns against a loop of single columns, with a small
   !!     budget that forces one sweep per column and one batch per chunk;
   !!  4. one column with T = 0 against xc_grad_potential_add, and the
   !!     symmetrisation of the input;
   !!  5. the refusals.
   !!
   !! Argument: the functional name (svwn, pbe, b3lyp).
   !!
   use pic_types, only: dp, int64
   use testdrive, only: new_unittest, unittest_type, error_type, check
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_linalg, only: trc_linalg_t
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_xc_functional, only: trc_xc_functional_t, xc_functional_by_name, XC_FAMILY_MGGA
   use trc_xc, only: trc_xc_rks, trc_xc_fxc_rks, trc_xc_result_t
   use trc_xc_grad, only: xc_grad_add, xc_grad_options, xc_grad_options_t
   use trc_xc_grad_potential, only: xc_grad_potential_add
   use trc_xc_grad_response, only: trc_xc_grad_response_many, trc_xc_grad_response_plan, xc_grad_response_core
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_dft_grid, only: dft_grid_t, dft_grid_build
   use trc_test_basis, only: read_xyz
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   implicit none
   private
   public :: collect_check_xc_grad_response, funcname_arg, release_cached

   character(len=*), parameter :: BASIS_JSON = "basis_sets/6-31gs.json"
   real(dp), parameter :: FD_H = 1.0e-3_dp
   !! Nuclear step of the Richardson pair (2h, h), Bohr.
   real(dp), parameter :: TOL_FD = 1.0e-9_dp
   !! Each column against the nuclear finite difference, absolute (Hartree/Bohr),
   !! over the atoms listed in fd_case (translation covers the rest).
   !! Measured worst 1.4e-10 (PBE, H2O2); the finite difference's own floor.
   real(dp), parameter :: TOL_ADD = 1.0e-10_dp
   !! The combined column against the sum of the two single-term columns.
   real(dp), parameter :: TOL_TRANS = 1.0e-9_dp
   !! The sum over atoms of a column.
   real(dp), parameter :: TOL_IDENT = 3.0e-7_dp
   !! The density-step identities, absolute. Measured worst 6.8e-8 (PBE): the t^4
   !! truncation left after the Richardson step.
   real(dp), parameter :: TOL_SAME = 1.0e-13_dp
   !! Batch against loop at the same chunking, relative to the largest element.
   !! The host measures the arithmetic per column as identical; a device build's
   !! atomics reorder sums, whose error stays well under this.
   real(dp), parameter :: TOL_SCREEN = 1.0e-9_dp
   !! What a 1e-12 density screen may drop from a column, absolute.
   real(dp), parameter :: SCREEN_COARSE = 1.0e-1_dp
   !! A screen coarse enough to drop batches that carry gradient.
   real(dp), parameter :: TOL_BATCH = 5.0e-13_dp
   !! Batch against loop, relative to the largest element. Measured worst 2.0e-13:
   !! a different column count changes the tile cap and so the grouping of
   !! the batches' partial sums, which is round-off over ~1e5 terms.
   real(dp), parameter :: TOL_POT = 1.0e-12_dp
   !! One column, T = 0, against xc_grad_potential_add, absolute.
   real(dp), parameter :: MIN_CONTROL = 1.0e-4_dp
   !! A negative control must miss by at least this, and a term be at least this large.

   character(len=16) :: funcname_arg = "svwn"

   type :: case_t
      !! One molecule's converged state, the reference grid and the probe densities.
      integer :: nat = 0, n = 0, nocc = 0
      integer, allocatable :: z(:)
      real(dp), allocatable :: r(:, :)
      type(trc_basis_t) :: bas
      type(trc_xc_functional_t) :: func
      type(trc_linalg_t) :: la
      type(dft_grid_t) :: grid
      type(trc_xc_grid_t) :: xg
      real(dp), allocatable :: dmat(:, :), cmo(:, :)
   end type case_t

   type(case_t), target, save :: cache(2)
   character(len=32), save :: cache_name(2) = ""

contains

   subroutine collect_check_xc_grad_response(testsuite)
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      testsuite = [ &
                  new_unittest("XC response gradient vs nuclear finite differences, water", test_fd_water), &
                  new_unittest("XC response gradient vs nuclear finite differences, H2O2", test_fd_h2o2), &
                  new_unittest("XC response gradient vs the density-step identities", test_identity), &
                  new_unittest("XC response gradient: batch == loop, small budget", test_batch), &
                  new_unittest("XC response gradient: no screen by default, screens change it", test_screen), &
                  new_unittest("XC response gradient: T = 0 is xc_grad_potential_add, symmetrised input", test_potential), &
                  new_unittest("XC response gradient refuses meta-GGA, VV10, range separation", test_refuse) &
                  ]
   end subroutine collect_check_xc_grad_response

   subroutine setup_case(cs, xyz, error)
      type(case_t), intent(inout) :: cs
      character(len=*), intent(in) :: xyz
      type(error_type), allocatable, intent(out) :: error
      type(trc_scf_options_t) :: opts
      type(trc_scf_result_t) :: res
      type(error_t) :: err
      integer :: dev
      real(dp), allocatable :: rr(:, :)

      dev = trc_bind_device(0)
      call read_xyz(trim(xyz), cs%nat, cs%z, rr)
      cs%r = rr
      call trc_basis_from_json(BASIS_JSON, cs%nat, cs%z, cs%r, cs%bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call cs%bas%to_device()
      cs%n = cs%bas%nao
      cs%nocc = nint(sum(real(cs%z, dp)))/2
      opts%conv_energy = 1.0e-12_dp
      opts%conv_density = 1.0e-11_dp
      opts%conv_diis = 1.0e-11_dp
      opts%max_iteration = 200
      opts%functional = trim(funcname_arg)
      opts%grid_level = 3
      call trc_scf_run(cs%bas, cs%nocc, cs%nocc, opts, res)
      call check(error, res%converged, "SCF did not converge")
      if (allocated(error)) return
      call check(error, xc_functional_by_name(trim(funcname_arg), cs%func), "unknown functional "//trim(funcname_arg))
      if (allocated(error)) return
      cs%dmat = res%dmat(:, :, 1)
      cs%cmo = res%cmo(:, :, 1)
      call dft_grid_build(cs%bas%at_r, cs%bas%at_elem, cs%grid, err, level=3)
      call check(error, .not. err%has_error(), "grid: "//err%get_message())
      if (allocated(error)) return
      call cs%xg%build_from_grid(cs%grid, cs%bas, max_pts=512)
      call cs%xg%to_device()
      call cs%la%init(cs%n)
   end subroutine setup_case

   subroutine acquire_case(cs, xyz, error)
      !! The converged state of `xyz`, built on first use and kept for the
      !! rest of the run (the SCF and the grid dominate a cheap check).
      type(case_t), pointer, intent(out) :: cs
      character(len=*), intent(in) :: xyz
      type(error_type), allocatable, intent(out) :: error
      integer :: k

      cs => null()
      do k = 1, size(cache)
         if (cache_name(k) == xyz) then
            cs => cache(k)
            return
         end if
      end do
      do k = 1, size(cache)
         if (len_trim(cache_name(k)) == 0) then
            call setup_case(cache(k), xyz, error)
            if (allocated(error)) return
            cache_name(k) = xyz
            cs => cache(k)
            return
         end if
      end do
      call check(error, .false., "case cache is full")
   end subroutine acquire_case

   subroutine release_cached()
      integer :: k
      do k = 1, size(cache)
         if (len_trim(cache_name(k)) > 0) call release_case(cache(k))
         cache_name(k) = ""
      end do
   end subroutine release_cached

   subroutine release_case(cs)
      type(case_t), intent(inout) :: cs
      call cs%la%release()
      call cs%xg%release()
      call cs%grid%destroy()
      call cs%bas%release()
   end subroutine release_case

   subroutine probe_densities(cs, p, pm, tm)
      !! Column p's relaxed-like density and transition density, in the MO span
      !! of the converged orbitals so their tails are the orbitals': P an
      !! occupied-occupied and virtual-virtual block of both signs (symmetric,
      !! indefinite, as a relaxed density is); T an occupied-virtual
      !! excitation, symmetrised. Phases of the orbitals do not matter: the
      !! same matrices are held fixed at every geometry.
      type(case_t), intent(in) :: cs
      integer, intent(in) :: p
      real(dp), intent(out) :: pm(cs%n, cs%n), tm(cs%n, cs%n)
      real(dp), allocatable :: co(:, :), cv(:, :), mo(:, :), mv(:, :), xov(:, :), tmp(:, :)
      integer :: i, j, nv, a, bq
      real(dp) :: ph

      nv = size(cs%cmo, 2) - cs%nocc
      ph = real(p, dp)
      co = cs%cmo(:, 1:cs%nocc)
      cv = cs%cmo(:, cs%nocc + 1:cs%nocc + nv)
      allocate (mo(cs%nocc, cs%nocc), mv(nv, nv), xov(cs%nocc, nv))
      do j = 1, cs%nocc
         do i = 1, cs%nocc
            mo(i, j) = 0.3_dp*cos(ph*real(i + 2*j, dp)) + 0.3_dp*cos(ph*real(j + 2*i, dp))
         end do
      end do
      do bq = 1, nv
         do a = 1, nv
            mv(a, bq) = 0.1_dp*cos(ph*real(a + 3*bq, dp)) + 0.1_dp*cos(ph*real(bq + 3*a, dp))
         end do
      end do
      do a = 1, nv
         do i = 1, cs%nocc
            xov(i, a) = 0.4_dp*cos(1.3_dp*ph*real(i, dp) + 0.7_dp*real(a, dp))
         end do
      end do
      pm = matmul(co, matmul(mo, transpose(co))) + matmul(cv, matmul(mv, transpose(cv)))
      tmp = matmul(co, matmul(xov, transpose(cv)))
      tm = tmp + transpose(tmp)
   end subroutine probe_densities

   subroutine analytic(cs, ncol, pm, tm, g, budget, nsweep, no_kxc, grid_response, screen_tol, tile_cap, error)
      !! xc_grad_response_core into a zeroed (3, nat, ncol) array.
      type(case_t), intent(inout) :: cs
      integer, intent(in) :: ncol
      real(dp), intent(in) :: pm(:, :, :), tm(:, :, :)
      real(dp), intent(out) :: g(3, cs%nat, ncol)
      integer(int64), intent(in), optional :: budget
      integer, intent(out), optional :: nsweep
      logical, intent(in), optional :: no_kxc, grid_response
      real(dp), intent(in), optional :: screen_tol
      integer(int64), intent(in), optional :: tile_cap
      type(error_type), allocatable, intent(out) :: error
      type(xc_grad_options_t) :: xopts
      type(error_t) :: err
      logical :: gr

      gr = .true.
      if (present(grid_response)) gr = grid_response
      xopts = xc_grad_options(grid_level=3, grid_max_pts=512, grid_response=gr, partition=cs%grid%partition)
      g = 0.0_dp
      !$acc enter data copyin(g)
      call xc_grad_response_core(cs%bas, cs%la, cs%func, cs%xg, cs%dmat, pm, tm, ncol, g, xopts, err, &
                                 budget=budget, nsweep=nsweep, screen_tol=screen_tol, no_kxc=no_kxc, tile_cap=tile_cap)
      !$acc update self(g)
      !$acc exit data delete(g)
      call check(error, .not. err%has_error(), err%get_message())
   end subroutine analytic

   subroutine terms_at(cs, rd, pm, tm, vals, error)
      !! (Tr(P V_xc[D]), 1/2 Tr(T f_xc[D](T))) with the basis and the grid
      !! rebuilt at the geometry rd (3, nat), the matrices held fixed.
      type(case_t), intent(in) :: cs
      real(dp), intent(in) :: rd(3, cs%nat)
      real(dp), intent(in) :: pm(cs%n, cs%n), tm(cs%n, cs%n)
      real(dp), intent(out) :: vals(2)
      type(error_type), allocatable, intent(out) :: error
      type(trc_basis_t) :: bas
      type(dft_grid_t) :: grid
      type(trc_xc_grid_t) :: xg
      type(trc_linalg_t) :: la
      type(trc_xc_result_t) :: xres
      type(error_t) :: err
      real(dp), allocatable :: vxc(:, :), vresp(:, :)

      allocate (vxc(cs%n, cs%n), vresp(cs%n, cs%n))
      call trc_basis_from_json(BASIS_JSON, cs%nat, cs%z, rd, bas, err, force_cart=.true.)
      call check(error, .not. err%has_error(), err%get_message())
      if (allocated(error)) return
      call bas%to_device()
      call dft_grid_build(bas%at_r, bas%at_elem, grid, err, level=3)
      call check(error, .not. err%has_error(), "grid: "//err%get_message())
      if (allocated(error)) return
      call xg%build_from_grid(grid, bas, max_pts=512)
      call xg%to_device()
      call la%init(cs%n)
      call trc_xc_rks(bas, xg, cs%func, cs%dmat, vxc, xres, la=la)
      vals(1) = sum(pm*vxc)
      call trc_xc_fxc_rks(bas, xg, cs%func, cs%dmat, tm, vresp, la=la)
      vals(2) = 0.5_dp*sum(tm*vresp)
      call la%release()
      call xg%release()
      call grid%destroy()
      call bas%release()
   end subroutine terms_at

   subroutine fd_gradient(cs, atoms, pm, tm, gfd, error)
      !! gfd(1:2, k, a): Richardson central differences of terms_at's two
      !! values with respect to coordinate k of atom a, for the atoms listed
      !! (the others stay zero: each atom costs 12 rebuilt grids).
      type(case_t), intent(in) :: cs
      integer, intent(in) :: atoms(:)
      real(dp), intent(in) :: pm(cs%n, cs%n), tm(cs%n, cs%n)
      real(dp), intent(out) :: gfd(2, 3, cs%nat)
      type(error_type), allocatable, intent(out) :: error
      real(dp) :: rd(3, cs%nat), vp(2), vm(2), vp2(2), vm2(2), d1(2), d2(2)
      integer :: a, k, ia

      gfd = 0.0_dp
      do ia = 1, size(atoms)
         a = atoms(ia)
         do k = 1, 3
            rd = cs%r
            rd(k, a) = cs%r(k, a) + FD_H
            call terms_at(cs, rd, pm, tm, vp, error)
            if (allocated(error)) return
            rd(k, a) = cs%r(k, a) - FD_H
            call terms_at(cs, rd, pm, tm, vm, error)
            if (allocated(error)) return
            rd(k, a) = cs%r(k, a) + 2.0_dp*FD_H
            call terms_at(cs, rd, pm, tm, vp2, error)
            if (allocated(error)) return
            rd(k, a) = cs%r(k, a) - 2.0_dp*FD_H
            call terms_at(cs, rd, pm, tm, vm2, error)
            if (allocated(error)) return
            d1 = (vp - vm)/(2.0_dp*FD_H)
            d2 = (vp2 - vm2)/(4.0_dp*FD_H)
            gfd(:, k, a) = (4.0_dp*d1 - d2)/3.0_dp
         end do
      end do
   end subroutine fd_gradient

   subroutine fd_case(xyz, label, atoms, error)
      !! The nuclear finite-difference check of one molecule: columns
      !! (P, 0), (0, T) and (P, T), each against its own finite difference.
      character(len=*), intent(in) :: xyz, label
      integer, intent(in) :: atoms(:)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), pointer :: cs
      real(dp), allocatable :: pm(:, :, :), tm(:, :, :), g(:, :, :), g_nk(:, :, :), g_ng(:, :, :), z(:, :)
      real(dp) :: gfd(2, 3, 4), dev(3), ref(3), trans, add, ctl_kxc, ctl_grid
      integer :: nc, c, a, k, ia
      character(len=240) :: line

      call acquire_case(cs, xyz, error)
      if (allocated(error)) return
      allocate (pm(cs%n, cs%n, 3), tm(cs%n, cs%n, 3), z(cs%n, cs%n))
      ! Column 1 is P alone, 2 is T alone, 3 both, from one probe.
      call probe_densities(cs, 1, z, tm(:, :, 2))
      call probe_densities(cs, 1, pm(:, :, 1), z)
      pm(:, :, 2) = 0.0_dp
      tm(:, :, 1) = 0.0_dp
      pm(:, :, 3) = pm(:, :, 1)
      tm(:, :, 3) = tm(:, :, 2)
      nc = 3
      allocate (g(3, cs%nat, nc), g_nk(3, cs%nat, nc), g_ng(3, cs%nat, nc))
      call analytic(cs, nc, pm, tm, g, error=error)
      if (allocated(error)) return
      call analytic(cs, nc, pm, tm, g_nk, no_kxc=.true., error=error)
      if (allocated(error)) return
      call analytic(cs, nc, pm, tm, g_ng, grid_response=.false., error=error)
      if (allocated(error)) return
      call fd_gradient(cs, atoms, pm(:, :, 3), tm(:, :, 3), gfd(:, :, 1:cs%nat), error)
      if (allocated(error)) return
      ! Columns 1 and 2 use the same P and T as column 3, so the three
      ! finite differences are the two terms and their sum.
      dev = 0.0_dp
      ref = 0.0_dp
      do c = 1, nc
         do ia = 1, size(atoms)
            a = atoms(ia)
            do k = 1, 3
               select case (c)
               case (1)
                  ref(1) = max(ref(1), abs(gfd(1, k, a)))
                  dev(1) = max(dev(1), abs(g(k, a, 1) - gfd(1, k, a)))
               case (2)
                  ref(2) = max(ref(2), abs(gfd(2, k, a)))
                  dev(2) = max(dev(2), abs(g(k, a, 2) - gfd(2, k, a)))
               case default
                  ref(3) = max(ref(3), abs(gfd(1, k, a) + gfd(2, k, a)))
                  dev(3) = max(dev(3), abs(g(k, a, 3) - gfd(1, k, a) - gfd(2, k, a)))
               end select
            end do
         end do
      end do
      add = maxval(abs(g(:, :, 3) - g(:, :, 1) - g(:, :, 2)))
      trans = maxval(abs(sum(g, dim=2)))
      ctl_kxc = 0.0_dp
      ctl_grid = 0.0_dp
      do ia = 1, size(atoms)
         a = atoms(ia)
         ctl_kxc = max(ctl_kxc, maxval(abs(g_nk(:, a, 2) - gfd(2, :, a))))
         ctl_grid = max(ctl_grid, maxval(abs(g_ng(:, a, 3) - gfd(1, :, a) - gfd(2, :, a))))
      end do
      write (line, "(a,a,a,3es10.2,a,3es10.2)") "  ", trim(label), ": max |analytic - FD| (P, T, P+T) ", dev, &
         "  max |term| ", ref
      print "(a)", trim(line)
      write (line, "(a,es10.2,a,es10.2,a,es10.2,a,es10.2)") "    additivity ", add, "  translation ", trans, &
         "  control without k_xc: ", ctl_kxc, "  control without grid response: ", ctl_grid
      print "(a)", trim(line)
      call check(error, all(ref > MIN_CONTROL), trim(label)//": a term is too small to test anything")
      if (.not. allocated(error)) call check(error, dev(1) < TOL_FD, trim(label)//": Tr(P V_xc) term differs from FD")
      if (.not. allocated(error)) call check(error, dev(2) < TOL_FD, trim(label)//": 1/2 f_xc T T term differs from FD")
      if (.not. allocated(error)) call check(error, dev(3) < TOL_FD, trim(label)//": combined column differs from FD")
      if (.not. allocated(error)) call check(error, add < TOL_ADD, trim(label)//": columns are not additive")
      if (.not. allocated(error)) call check(error, trans < TOL_TRANS, trim(label)//": not translationally invariant")
      if (.not. allocated(error)) call check(error, ctl_kxc > MIN_CONTROL, &
                                             trim(label)//": dropping k_xc does not break the T term")
      if (.not. allocated(error)) call check(error, ctl_grid > 100.0_dp*TOL_FD, &
                                             trim(label)//": dropping the grid response does not break the gradient")
   end subroutine fd_case

   subroutine test_fd_water(error)
      type(error_type), allocatable, intent(out) :: error
      call fd_case("water.xyz", "water", [2], error)
   end subroutine test_fd_water

   subroutine test_fd_h2o2(error)
      type(error_type), allocatable, intent(out) :: error
      call fd_case("h2o2.xyz", "h2o2", [1], error)
   end subroutine test_fd_h2o2

   subroutine xc_gradient_at(cs, dens, grad)
      !! trc_xc_grad's own E_xc gradient (grid response on) at the density
      !! `dens`, on its own level-3 grid.
      type(case_t), intent(inout) :: cs
      real(dp), intent(in) :: dens(:, :)
      real(dp), intent(out) :: grad(:, :)
      real(dp), allocatable :: d3(:, :, :)
      allocate (d3(size(dens, 1), size(dens, 2), 1))
      d3(:, :, 1) = dens
      grad = 0.0_dp
      !$acc enter data copyin(d3, grad)
      call xc_grad_add(cs%bas, cs%la, cs%func, 1, d3, grad, xc_grad_options(grid_level=3, grid_max_pts=512))
      !$acc update self(grad)
      !$acc exit data delete(d3, grad)
   end subroutine xc_gradient_at

   subroutine test_identity(error)
      !! d/dR of the T term is 1/2 d^2/dt^2 of trc_xc_grad's gradient at
      !! D + t T (E_xc[D + tT] = E + t E' + t^2/2 E'' + ...), and of the P term
      !! d/dt, each Richardson-extrapolated over t and t/2. No nuclear
      !! difference: an independent route to the same derivative through the
      !! existing, separately validated ground-state gradient.
      type(error_type), allocatable, intent(out) :: error
      type(case_t), pointer :: cs
      real(dp), allocatable :: pm(:, :, :), tm(:, :, :), g(:, :, :)
      real(dp), allocatable :: g0(:, :), gp(:, :), gm(:, :), gp2(:, :), gm2(:, :), s1(:, :), s2(:, :)
      real(dp) :: t, dev_t, dev_p, mag_t
      integer :: nat

      call acquire_case(cs, "water.xyz", error)
      if (allocated(error)) return
      nat = cs%nat
      allocate (pm(cs%n, cs%n, 2), tm(cs%n, cs%n, 2), g(3, nat, 2))
      allocate (g0(3, nat), gp(3, nat), gm(3, nat), gp2(3, nat), gm2(3, nat), s1(3, nat), s2(3, nat))
      pm = 0.0_dp
      tm = 0.0_dp
      call probe_densities(cs, 1, pm(:, :, 1), tm(:, :, 2))
      ! Column 1: P alone; column 2: T alone.
      tm(:, :, 1) = 0.0_dp
      call analytic(cs, 2, pm, tm, g, error=error)
      if (allocated(error)) return
      ! tm(:,:,2) came from the same call that filled pm(:,:,2) with the P: only T is kept there.
      t = 0.05_dp
      call xc_gradient_at(cs, cs%dmat, g0)
      ! The T term, second difference.
      call xc_gradient_at(cs, cs%dmat + t*tm(:, :, 2), gp)
      call xc_gradient_at(cs, cs%dmat - t*tm(:, :, 2), gm)
      call xc_gradient_at(cs, cs%dmat + 0.5_dp*t*tm(:, :, 2), gp2)
      call xc_gradient_at(cs, cs%dmat - 0.5_dp*t*tm(:, :, 2), gm2)
      s1 = (gp + gm - 2.0_dp*g0)/(t*t)
      s2 = (gp2 + gm2 - 2.0_dp*g0)/(0.25_dp*t*t)
      dev_t = maxval(abs(g(:, :, 2) - 0.5_dp*(4.0_dp*s2 - s1)/3.0_dp))
      mag_t = maxval(abs(g(:, :, 2)))
      ! The P term, first difference, at a smaller step: P has a
      ! virtual-virtual block, so D + tP leaves the density's tails sooner.
      t = 0.01_dp
      call xc_gradient_at(cs, cs%dmat + t*pm(:, :, 1), gp)
      call xc_gradient_at(cs, cs%dmat - t*pm(:, :, 1), gm)
      call xc_gradient_at(cs, cs%dmat + 0.5_dp*t*pm(:, :, 1), gp2)
      call xc_gradient_at(cs, cs%dmat - 0.5_dp*t*pm(:, :, 1), gm2)
      s1 = (gp - gm)/(2.0_dp*t)
      s2 = (gp2 - gm2)/t
      dev_p = maxval(abs(g(:, :, 1) - (4.0_dp*s2 - s1)/3.0_dp))
      print "(a,es10.2,a,es10.2,a,es10.2)", "  identity: T term max |analytic - 1/2 d2/dt2| ", dev_t, &
         "  P term ", dev_p, "  max |T term| ", mag_t
      call check(error, mag_t > MIN_CONTROL, "the T term is too small to test anything")
      if (.not. allocated(error)) call check(error, dev_t < TOL_IDENT, "T term differs from the second density difference")
      if (.not. allocated(error)) call check(error, dev_p < TOL_IDENT, "P term differs from the density difference")
   end subroutine test_identity

   subroutine test_batch(error)
      !! Batch against loop, three ways. Same chunking (one tile cap for both
      !! sides): only the per-column arithmetic differs, so near bitwise on
      !! the host. Different chunking (default plan, a one-column budget, a
      !! two-column budget with several chunks): the partial sums regroup,
      !! round-off at the scaled tolerance.
      type(error_type), allocatable, intent(out) :: error
      type(case_t), pointer :: cs
      real(dp), allocatable :: pm(:, :, :), tm(:, :, :), gb(:, :, :), gl(:, :, :), g1(:, :, :), gs(:, :, :), gt(:, :, :)
      real(dp), allocatable :: gq(:, :, :), gc(:, :, :)
      integer(int64) :: tcap_small, tcap_big, tcap2, need, bud, bfound
      integer :: nc_small, nc_big, nsw_small, nsw_big, c, nsw1, nc2, nch2, nch_big, nsw2, it, nsw_same
      real(dp) :: scale, dev_loop, dev_budget, dev_same, dev_two
      integer, parameter :: NCOL = 4

      call acquire_case(cs, "water.xyz", error)
      if (allocated(error)) return
      allocate (pm(cs%n, cs%n, NCOL), tm(cs%n, cs%n, NCOL))
      allocate (gb(3, cs%nat, NCOL), gl(3, cs%nat, NCOL), g1(3, cs%nat, 1), gs(3, cs%nat, NCOL), gt(3, cs%nat, NCOL))
      allocate (gq(3, cs%nat, 3), gc(3, cs%nat, NCOL))
      do c = 1, NCOL
         call probe_densities(cs, c, pm(:, :, c), tm(:, :, c))
      end do
      call trc_xc_grad_response_plan(cs%bas, cs%xg, NCOL, nc=nc_big, tcap=tcap_big, nchunk=nch_big)
      call analytic(cs, NCOL, pm, tm, gb, nsweep=nsw_big, error=error)
      if (allocated(error)) return
      ! A budget of one doubles: one sweep per column, the tile cap at its floor
      ! of one batch per chunk.
      call trc_xc_grad_response_plan(cs%bas, cs%xg, NCOL, 1_int64, nc=nc_small, tcap=tcap_small)
      call analytic(cs, NCOL, pm, tm, gs, budget=1_int64, nsweep=nsw_small, error=error)
      if (allocated(error)) return
      ! Same chunking: the small tile cap on both sides.
      call analytic(cs, NCOL, pm, tm, gt, tile_cap=tcap_small, nsweep=nsw_same, error=error)
      if (allocated(error)) return
      do c = 1, NCOL
         call analytic(cs, 1, pm(:, :, c:c), tm(:, :, c:c), g1, tile_cap=tcap_small, nsweep=nsw1, error=error)
         if (allocated(error)) return
         gc(:, :, c) = g1(:, :, 1)
      end do
      dev_same = maxval(abs(gt - gc))/maxval(abs(gc))
      ! Default chunking on both sides (the tile cap follows the column count).
      do c = 1, NCOL
         call analytic(cs, 1, pm(:, :, c:c), tm(:, :, c:c), g1, nsweep=nsw1, error=error)
         if (allocated(error)) return
         gl(:, :, c) = g1(:, :, 1)
      end do
      scale = maxval(abs(gl))
      dev_loop = maxval(abs(gb - gl))/scale
      dev_budget = maxval(abs(gs - gl))/scale
      ! Two columns per sweep and several chunks per sweep: the budget is
      ! searched downward from the three-column plan's need until the plan
      ! says nc = 2 with more than one chunk.
      call trc_xc_grad_response_plan(cs%bas, cs%xg, 3, nc=nc2, tcap=tcap2, nchunk=nch2, need=need)
      bfound = 0
      bud = need
      do it = 1, 200
         bud = (bud*19_int64)/20_int64
         call trc_xc_grad_response_plan(cs%bas, cs%xg, 3, bud, nc=nc2, tcap=tcap2, nchunk=nch2)
         if (nc2 == 2 .and. nch2 > 1) then
            bfound = bud
            exit
         end if
         if (nc2 < 2) exit
      end do
      dev_two = -1.0_dp
      nsw2 = 0
      if (bfound > 0) then
         call analytic(cs, 3, pm(:, :, 1:3), tm(:, :, 1:3), gq, budget=bfound, nsweep=nsw2, error=error)
         if (allocated(error)) return
         dev_two = maxval(abs(gq - gl(:, :, 1:3)))/scale
      end if
      print "(a,es10.2,a,es10.2,a,es10.2,a,es10.2)", "  batch: rel |same chunking| ", dev_same, "  |one sweep - loop| ", &
         dev_loop, "  |one-column budget - loop| ", dev_budget, "  |nc=2, several chunks - loop| ", dev_two
      print "(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)", "    sweeps ", nsw_big, " and ", nsw_small, "; same-chunking sweeps ", nsw_same, &
         "; chunks/sweep at default ", nch_big, "; nc=2 case: chunks ", nch2, " sweeps ", nsw2
      call check(error, scale > MIN_CONTROL, "the gradients are too small to test anything")
      if (.not. allocated(error)) call check(error, nsw_big == 1 .and. nc_big == NCOL, "the default budget must take all columns in one sweep")
      if (.not. allocated(error)) call check(error, nsw_small == NCOL .and. nc_small == 1, &
                                             "a one-column budget must take one sweep per column")
      if (.not. allocated(error)) call check(error, tcap_small < tcap_big, "the small budget must shrink the tile cap")
      if (.not. allocated(error)) call check(error, nsw_same == 1, "the same-chunking batch must be one sweep")
      if (.not. allocated(error)) call check(error, bfound > 0 .and. nsw2 == 2, &
                                             "no budget gives two columns per sweep with several chunks")
      if (.not. allocated(error)) call check(error, dev_same < TOL_SAME, "same chunking: batch differs from the loop")
      if (.not. allocated(error)) call check(error, dev_loop < TOL_BATCH, "one sweep differs from the loop of single columns")
      if (.not. allocated(error)) call check(error, dev_budget < TOL_BATCH, "a small budget changes the answer")
      if (.not. allocated(error)) call check(error, dev_two < TOL_BATCH, "two columns per sweep with several chunks change the answer")
   end subroutine test_batch

   subroutine test_screen(error)
      !! The default is no screen, the operator trc_xc_fxc_rks and
      !! trc_xc_kxc_rks_apply define. A 1e-12 screen (xc_grad_potential_add's)
      !! drops a negligible part; a coarse one visibly changes the gradient,
      !! so the choice is exercised: the default agrees with the unscreened
      !! finite differences of fd_case, the coarse screen would not.
      type(error_type), allocatable, intent(out) :: error
      type(case_t), pointer :: cs
      real(dp), allocatable :: pm(:, :, :), tm(:, :, :), g0(:, :, :), g12(:, :, :), g4(:, :, :)
      real(dp) :: d12, d4

      call acquire_case(cs, "water.xyz", error)
      if (allocated(error)) return
      allocate (pm(cs%n, cs%n, 1), tm(cs%n, cs%n, 1), g0(3, cs%nat, 1), g12(3, cs%nat, 1), g4(3, cs%nat, 1))
      call probe_densities(cs, 1, pm(:, :, 1), tm(:, :, 1))
      call analytic(cs, 1, pm, tm, g0, error=error)
      if (allocated(error)) return
      call analytic(cs, 1, pm, tm, g12, screen_tol=1.0e-12_dp, error=error)
      if (allocated(error)) return
      call analytic(cs, 1, pm, tm, g4, screen_tol=SCREEN_COARSE, error=error)
      if (allocated(error)) return
      d12 = maxval(abs(g0 - g12))
      d4 = maxval(abs(g0 - g4))
      print "(a,es10.2,a,es10.2)", "  screen: max |unscreened - 1e-12 screen| ", d12, "  |unscreened - coarse screen| ", d4
      call check(error, d12 < TOL_SCREEN, "the 1e-12 screen drops more than it should")
      if (.not. allocated(error)) call check(error, d4 > 100.0_dp*TOL_FD, "a coarse screen does not change the gradient")
   end subroutine test_screen

   subroutine test_potential(error)
      type(error_type), allocatable, intent(out) :: error
      type(case_t), pointer :: cs
      type(xc_grad_options_t) :: xopts
      real(dp), allocatable :: pm(:, :, :), tm(:, :, :), g(:, :, :), gref(:, :), pns(:, :, :), tns(:, :, :), gn(:, :, :)
      real(dp) :: dev, dev_sym
      integer :: n

      call acquire_case(cs, "water.xyz", error)
      if (allocated(error)) return
      n = cs%n
      allocate (pm(n, n, 1), tm(n, n, 1), g(3, cs%nat, 1), gref(3, cs%nat), gn(3, cs%nat, 1))
      call probe_densities(cs, 3, pm(:, :, 1), tm(:, :, 1))
      tm = 0.0_dp
      call analytic(cs, 1, pm, tm, g, screen_tol=1.0e-12_dp, error=error)
      if (allocated(error)) return
      xopts = xc_grad_options(grid_level=3, grid_max_pts=512, grid_response=.true., partition=cs%grid%partition)
      gref = 0.0_dp
      !$acc enter data copyin(gref)
      call xc_grad_potential_add(cs%bas, cs%la, cs%func, cs%xg, cs%dmat, pm(:, :, 1), gref, xopts)
      !$acc update self(gref)
      !$acc exit data delete(gref)
      dev = maxval(abs(g(:, :, 1) - gref))
      ! Non-symmetric input: P + A and T + A' with A antisymmetric give the same answer.
      allocate (pns(n, n, 1), tns(n, n, 1))
      call probe_densities(cs, 2, pm(:, :, 1), tm(:, :, 1))
      call skew(pm(:, :, 1), pns(:, :, 1))
      call skew(tm(:, :, 1), tns(:, :, 1))
      call analytic(cs, 1, pm, tm, g, error=error)
      if (allocated(error)) return
      call analytic(cs, 1, pns, tns, gn, error=error)
      if (allocated(error)) return
      dev_sym = maxval(abs(g - gn))
      print "(a,es10.2,a,es10.2,a,es10.2)", "  T = 0 vs xc_grad_potential_add: max |diff| ", dev, &
         "  unsymmetrised input: max |diff| ", dev_sym, "  max |term| ", maxval(abs(gref))
      call check(error, maxval(abs(gref)) > MIN_CONTROL, "the term is too small to test anything")
      if (.not. allocated(error)) call check(error, dev < TOL_POT, "T = 0 differs from xc_grad_potential_add")
      if (.not. allocated(error)) call check(error, dev_sym < TOL_POT, "an unsymmetrised input changes the answer")
   end subroutine test_potential

   subroutine skew(x, y)
      !! y = x + a deterministic antisymmetric matrix (x symmetric): the same symmetric part.
      real(dp), intent(in) :: x(:, :)
      real(dp), intent(out) :: y(:, :)
      integer :: i, j
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            y(i, j) = x(i, j) + 0.1_dp*sin(0.3_dp*real(i*j, dp))*sign(1.0_dp, real(i - j, dp))
            if (i == j) y(i, j) = x(i, j)
         end do
      end do
   end subroutine skew

   subroutine test_refuse(error)
      !! No SCF needed: the refusal precedes any use of the basis or grid.
      type(error_type), allocatable, intent(out) :: error
      type(trc_xc_functional_t) :: f_mgga, f_rsh, f_vv10
      type(trc_basis_t) :: bas
      type(trc_xc_grid_t) :: xg
      type(trc_linalg_t) :: la
      type(xc_grad_options_t) :: xopts
      type(error_t) :: err
      real(dp) :: d(0, 0), p(0, 0, 1), t(0, 0, 1), g(3, 0, 1)

      f_mgga%family = XC_FAMILY_MGGA
      f_rsh%omega = 0.3_dp
      f_rsh%exx_lr = 0.8_dp
      f_vv10%nlc_b = 6.0_dp
      f_vv10%nlc_c = 0.01_dp
      call trc_xc_grad_response_many(bas, la, f_mgga, xg, d, p, t, 1, g, xopts, err)
      call check(error, err%has_error() .and. index(err%get_message(), "meta-GGA") > 0, "meta-GGA must be refused")
      if (allocated(error)) return
      call trc_xc_grad_response_many(bas, la, f_rsh, xg, d, p, t, 1, g, xopts, err)
      call check(error, err%has_error() .and. index(err%get_message(), "range-separated") > 0, "RSH must be refused")
      if (allocated(error)) return
      call trc_xc_grad_response_many(bas, la, f_vv10, xg, d, p, t, 1, g, xopts, err)
      call check(error, err%has_error() .and. index(err%get_message(), "VV10") > 0, "VV10 must be refused")
   end subroutine test_refuse

end module check_xc_grad_response_suite

program check_xc_grad_response
   use, intrinsic :: iso_fortran_env, only: error_unit
   use testdrive, only: run_testsuite
   use check_xc_grad_response_suite, only: collect_check_xc_grad_response, funcname_arg, release_cached
   implicit none
   integer :: stat

   if (command_argument_count() >= 1) call get_command_argument(1, funcname_arg)
   stat = 0
   call run_testsuite(collect_check_xc_grad_response, error_unit, stat, parallel=.false.)
   call release_cached()
   if (stat > 0) then
      write (error_unit, "(i0, a)") stat, " test(s) failed"
      error stop 1
   end if
end program check_xc_grad_response
