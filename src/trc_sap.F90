!! The superposition-of-atomic-potentials guess
module trc_sap
   !! SAP builds a one-electron Hamiltonian out of free-atom potentials and
   !! diagonalises it, rather than superposing free-atom densities and
   !! diagonalising the Fock matrix that density makes. Lehtola's point (JCTC
   !! 15 1593, 2019) is that the atomic POTENTIAL is the better thing to
   !! superpose: it is what an electron in the molecule actually sees, and it
   !! carries the atomic shell structure into the guess orbitals directly
   !! instead of through one Fock build.
   !!
   !!   H_sap = T + V_ne + sum_A Vscr_A(|r - R_A|)
   !!
   !! TWO HONEST QUALIFICATIONS.
   !!
   !! First, this is not the cheap SAP. Lehtola's is cheap because the atomic
   !! potentials are TABULATED, so the guess needs no atomic calculation at
   !! all. There is no such table in this tree and inventing one is not a
   !! thing to do quietly, so `Vscr_A` is derived from the same free-atom SCF
   !! that SAD already runs and caches. SAP therefore costs SAD plus a
   !! quadrature: it is here for guess QUALITY, not for speed. Dropping a
   !! real table in later would make it cheaper than SAD, and `atom_screen`
   !! is the only routine that would have to change.
   !!
   !! Second, the atomic effective potential of a Hartree-Fock atom is
   !! nonlocal, and this uses a local approximation to it -- the atom's own
   !! Hartree potential plus the Dirac-Slater exchange potential of its
   !! density. That is the same approximation the tabulated potentials embody
   !! and it is the reason SAP is a guess and not a method.
   !!
   !! WHAT IS ON THE GRID, AND WHAT IS NOT. The nuclear attraction -Z/r is
   !! singular and is already available exactly, so it is NOT quadratured: T
   !! and V_ne come in as the analytic one-electron matrices and only the
   !! SCREENING potential, which is smooth and short ranged, is integrated.
   !! Putting -Z/r on a grid instead would be both slower and worse.
   use trc_boys, only: dp
   use trc_api, only: trc_basis_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_sad, only: trc_sad_build, trc_sac_build, trc_atom_ranges
   use trc_dft_radial, only: treutler_ahlrichs_radial
   use trc_lebedev, only: lebedev_grid, lebedev_order_at_least
   use trc_dft_grid, only: dft_grid_t, build_dft_grid
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_collocation, only: shell_collocate, NCART_MAX
   use pic_blas_interfaces, only: pic_gemm
   implicit none
   private

   public :: trc_sap_build, trc_sap_screen_table

   !> Radial nodes per element for the screening table. The potential is
   !> smooth, and the M4 mapping puts the nodes where it bends, so this is
   !> generous rather than tuned; the table is built once per element.
   integer, parameter :: SAP_NRAD = 400

   !> Angular order for the spherical average of the atomic density. The
   !> free-atom density is already spherical -- fractional occupations over
   !> the degenerate frontier are what make it so -- but averaging costs
   !> nothing here and means a non-spherical atom degrades the guess instead
   !> of silently tilting it.
   integer, parameter :: SAP_NANG = 26

   real(dp), parameter :: PI = 3.14159265358979323846_dp

   !> Points per quadrature batch for the assembly.
   integer, parameter :: SAP_BATCH = 128

   !
   !> Quadrature level for the screening integral.
   !>
   !> This is the DFT default, and it is NOT a lazy inheritance -- a cheaper
   !> grid was tried first, on the reasoning that a guess thrown away after
   !> one diagonalisation cannot need a grid chosen to converge an
   !> exchange-correlation energy. That reasoning is wrong. Quadrature error
   !> in a matrix element is not a small perturbation of the guess, it is a
   !> wrong operator, and the orbitals it produces are worse rather than
   !> approximate. On the 123-atom silica slice, cc-pVDZ, at iteration 12:
   !>
   !>              E(12)                 dE(12)     guess cost
   !>   level 0    -11194.245761258     -9.5e-05      7 s
   !>   level 1    -11194.245769153     -9.9e-05     25 s
   !>   level 3    -11194.245774662     -3.3e-10     80 s
   !>   SAD        -11194.245774661     -6.0e-09      <1 s
   !>
   !> Level 3 is a better guess than SAD, by about one iteration. Levels 0
   !> and 1 are not converging at all by iteration 12 -- worse than SAD and
   !> cheaper than nothing is worth.
   !>
   !> WHICH MEANS SAP IS CURRENTLY A LOSS AT THIS SIZE. The one iteration it
   !> saves is ~12 s and the guess costs 80, of which 18 is `assemble` and 4
   !> the per-point potential, both host-side and single-threaded. The fix is
   !> to put those on the device the way trc_xc does, not to coarsen the
   !> grid; until then SAP earns its place on small molecules (water/6-31G:
   !> 9 iterations against 10, at no measurable cost) and not on slabs.
   !>
   !> TRC_SAP_LEVEL overrides it so the trade can be re-measured without a
   !> rebuild, and an explicit `level` argument overrides both.
   !
   integer, parameter :: SAP_GRID_LEVEL = 3

contains

   pure integer function ncart(l)
      integer, intent(in) :: l
      ncart = (l + 1)*(l + 2)/2
   end function ncart

   !
   ! The radial screening potential of one free atom, on that element's own
   ! radial mesh.
   !
   !   Vscr(r) = V_H[rho](r) + v_x[rho](r)
   !
   ! with rho the spherically averaged free-atom density and
   !
   !   V_H(r) = 4 pi [ (1/r) int_0^r s^2 rho ds + int_r^inf s rho ds ]
   !   v_x(r) = -(3 rho / pi)^(1/3)                      Dirac-Slater
   !
   ! `nel` comes back as 4 pi int r^2 rho dr, the electron count the table
   ! actually integrates to. It is both a check on the quadrature and the
   ! coefficient of the r -> infinity tail: past the last node the atom looks
   ! like a point charge of -nel, so Vscr = nel/r, and V_ne + Vscr then goes
   ! to zero for a neutral atom instead of to a constant.
   !
   subroutine trc_sap_screen_table(b, ish0, ish1, centre, z, dblk, r, vscr, nel, error)
      type(trc_basis_t), intent(in) :: b
      integer,  intent(in) :: ish0, ish1     !! the atom's shell range in `b`
      real(dp), intent(in) :: centre(3)
      integer,  intent(in) :: z
      real(dp), intent(in) :: dblk(:, :)     !! its density, in those functions
      real(dp), allocatable, intent(out) :: r(:), vscr(:)
      real(dp), intent(out) :: nel
      type(error_t), intent(inout) :: error

      real(dp), allocatable :: dr(:), rho(:), ang(:, :), wang(:)
      real(dp) :: chi(NCART_MAX), gchi(3, NCART_MAX)
      real(dp), allocatable :: phi(:)
      real(dp) :: pt(3), d(3), s, acc_in, acc_out, rh
      integer :: nrad, nang, k, ia, ish, np, nc, off, nz, i, j

      nz = size(dblk, 1)
      call treutler_ahlrichs_radial(SAP_NRAD, z, r, dr, error)
      if (error%has_error()) return
      nrad = size(r)
      nang = lebedev_order_at_least(SAP_NANG)
      call lebedev_grid(nang, ang, wang, error)
      if (error%has_error()) return

      allocate (rho(nrad), phi(nz), vscr(nrad))

      ! --- the spherically averaged density on the radial mesh -------------
      do k = 1, nrad
         rh = 0.0_dp
         do ia = 1, nang
            pt = centre + r(k)*ang(:, ia)
            phi = 0.0_dp
            do ish = ish0, ish1
               np = b%sh_np(ish)
               nc = ncart(b%sh_l(ish))
               d = pt - b%sh_r(:, ish)
               call shell_collocate(b%sh_l(ish), np, b%sh_e(1:np, ish), b%sh_c(1:np, ish), d, chi, gchi)
               off = b%sh_ao(ish) - b%sh_ao(ish0)
               do i = 1, nc
                  phi(off + i) = chi(i)
               end do
            end do
            ! rho = phi^T D phi, and wang sums to one, so this is the average
            s = 0.0_dp
            do j = 1, nz
               do i = 1, nz
                  s = s + phi(i)*dblk(i, j)*phi(j)
               end do
            end do
            rh = rh + wang(ia)*s
         end do
         rho(k) = max(rh, 0.0_dp)
      end do

      ! --- the two radial integrals, and the electron count ----------------
      ! acc_in(k) = int_0^r_k s^2 rho ds is accumulated forward; the outward
      ! one has to be accumulated backward, so it is done in a second sweep
      ! rather than by subtracting, which would lose the small tail to
      ! cancellation against the large head.
      nel = 0.0_dp
      do k = 1, nrad
         nel = nel + r(k)*r(k)*rho(k)*dr(k)
      end do
      nel = 4.0_dp*PI*nel

      acc_in = 0.0_dp
      do k = 1, nrad
         acc_in = acc_in + r(k)*r(k)*rho(k)*dr(k)
         vscr(k) = acc_in/r(k)
      end do
      acc_out = 0.0_dp
      do k = nrad, 1, -1
         vscr(k) = 4.0_dp*PI*(vscr(k) + acc_out)
         acc_out = acc_out + r(k)*rho(k)*dr(k)
      end do

      ! --- Dirac-Slater exchange ------------------------------------------
      do k = 1, nrad
         if (rho(k) > 0.0_dp) vscr(k) = vscr(k) - (3.0_dp*rho(k)/PI)**(1.0_dp/3.0_dp)
      end do
      deallocate (dr, rho, ang, wang, phi)
   end subroutine trc_sap_screen_table

   !
   ! Vscr of one element at an arbitrary distance.
   !
   ! Linear between nodes; nel/rr past the last one, where the atom is a
   ! point charge; the first node's value inside it, which is a region no
   ! molecular grid resolves and where the nuclear term dominates anyway.
   !
   pure real(dp) function interp(r, v, nel, rr)
      real(dp), intent(in) :: r(:), v(:), nel, rr
      integer :: lo, hi, mid
      real(dp) :: t
      if (rr <= r(1)) then
         interp = v(1)
      else if (rr >= r(size(r))) then
         interp = nel/rr
      else
         lo = 1; hi = size(r)
         do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (r(mid) > rr) then
               hi = mid
            else
               lo = mid
            end if
         end do
         t = (rr - r(lo))/(r(hi) - r(lo))
         interp = v(lo)*(1.0_dp - t) + v(hi)*t
      end if
   end function interp

   !
   ! H_sap = T + V_ne + sum_A Vscr_A.
   !
   ! `tmat` and `vmat` are the analytic one-electron matrices the caller
   ! already has from trc_1e. `nelec` present selects SAC atoms -- charge
   ! spread over the atoms -- rather than SAD's neutral ones, which is worth
   ! having for an ion for the same reason it is worth having in SAC.
   !
   subroutine trc_sap_build(b, tmat, vmat, hsap, error, verbose, nelec, level)
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: tmat(b%nao, b%nao), vmat(b%nao, b%nao)
      real(dp), intent(out) :: hsap(b%nao, b%nao)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose
      integer, intent(in), optional :: nelec   !! molecular electron count, for SAC atoms
      integer, intent(in), optional :: level   !! quadrature level, as build_dft_grid spells it

      real(dp), allocatable :: dguess(:, :)
      integer,  allocatable :: atom_of(:), first(:), last(:), zint(:)
      type dtab
         real(dp), allocatable :: r(:), v(:)
         real(dp) :: nel = 0.0_dp
         integer :: z = 0, nz = 0
      end type dtab
      type(dtab), allocatable :: tab(:)
      integer, allocatable :: tab_of(:)
      type(dft_grid_t) :: g
      real(dp), allocatable :: vpt(:)
      real(dp) :: rr
      integer :: ia, ib, k, a0, a1, nz, ntab, ip, lv
      integer :: c0, c1, crate
      logical :: talk

      talk = .false.
      if (present(verbose)) talk = verbose
      hsap = 0.0_dp
      call system_clock(c0, crate)

      ! --- the free atoms, and their ranges -------------------------------
      if (present(nelec)) then
         call trc_sac_build(b, nelec, dguess, error, verbose=verbose)
      else
         call trc_sad_build(b, dguess, error, verbose=verbose)
      end if
      if (error%has_error()) return
      call trc_atom_ranges(b, atom_of, first, last, error)
      if (error%has_error()) return
      call stage('free atoms   ')

      ! --- one screening table per distinct (element, function count) ------
      allocate (tab(b%natm), tab_of(b%natm), zint(b%natm))
      tab_of = 0
      ntab = 0
      do ia = 1, b%natm
         if (first(ia) == 0) cycle
         zint(ia) = nint(b%at_z(ia))
         a0 = b%sh_ao(first(ia))
         a1 = b%sh_ao(last(ia)) + ncart(b%sh_l(last(ia))) - 1
         nz = a1 - a0 + 1
         do k = 1, ntab
            if (tab(k)%z == zint(ia) .and. tab(k)%nz == nz) tab_of(ia) = k
         end do
         if (tab_of(ia) == 0) then
            ntab = ntab + 1
            call trc_sap_screen_table(b, first(ia), last(ia), b%at_r(:, ia), zint(ia), &
                                      dguess(a0:a1, a0:a1), tab(ntab)%r, tab(ntab)%v, &
                                      tab(ntab)%nel, error)
            if (error%has_error()) return
            tab(ntab)%z = zint(ia); tab(ntab)%nz = nz
            if (talk) print '(a,i0,a,f12.8,a,f12.8)', "  sap: Z = ", zint(ia), &
               "  radial density integrates to ", tab(ntab)%nel, &
               "   r*Vscr at the last node ", tab(ntab)%r(size(tab(ntab)%r))* &
               tab(ntab)%v(size(tab(ntab)%v))
            tab_of(ia) = ntab
         end if
      end do

      call stage('screen tables')

      ! --- the molecular quadrature ---------------------------------------
      lv = SAP_GRID_LEVEL
      if (present(level)) lv = level
      block
         character(len=16) :: e
         integer :: ios, lv_env
         e = ' '
         call get_environment_variable('TRC_SAP_LEVEL', e)
         if (len_trim(e) > 0 .and. .not. present(level)) then
            read (e, *, iostat=ios) lv_env
            if (ios == 0) lv = lv_env
         end if
      end block
      call build_dft_grid(b%at_r(:, 1:b%natm), zint, g, error, level=lv)
      if (error%has_error()) return
      if (talk) print '(a,i0,a,i0,a)', "  sap: quadrature level ", lv, ", ", &
         g%n_points, " points"

      ! The potential is a scalar on the grid, and every atom contributes to
      ! every point: nothing is screened here because Vscr falls off as
      ! nel/r, which is not small at molecular distances -- it is the
      ! long-range neutrality of V_ne + Vscr, and dropping it would leave a
      ! spurious charge behind.
      call stage('grid         ')
      allocate (vpt(g%n_points))
      do ip = 1, g%n_points
         vpt(ip) = 0.0_dp
         do ia = 1, b%natm
            if (tab_of(ia) == 0) cycle
            k = tab_of(ia)
            rr = sqrt(sum((g%coords(:, ip) - b%at_r(:, ia))**2))
            vpt(ip) = vpt(ip) + interp(tab(k)%r, tab(k)%v, tab(k)%nel, rr)
         end do
         vpt(ip) = vpt(ip)*g%weights(ip)
      end do

      call stage('potential    ')
      call assemble(b, g, vpt, hsap, error)
      if (error%has_error()) return
      call stage('assemble     ')

      hsap = hsap + tmat + vmat
      call g%destroy()
      deallocate (dguess, atom_of, first, last, zint, tab, tab_of, vpt)

   contains

      !> Stage timing, printed only when TRC_BUILD_TIMING is set, as the Fock
      !> build's own stages are. Four routines with very different costs and
      !> no way to tell from the outside which one is slow.
      subroutine stage(what)
         character(len=*), intent(in) :: what
         character(len=8) :: e
         e = ' '
         call get_environment_variable('TRC_BUILD_TIMING', e)
         call system_clock(c1)
         if (len_trim(e) > 0) print '(a,a,f9.3,a)', '  [sap] ', what, &
            real(c1 - c0, dp)/real(crate, dp), ' s'
         c0 = c1
      end subroutine stage
   end subroutine trc_sap_build

   !
   ! H_munu += sum_g wv(g) chi_mu(g) chi_nu(g).
   !
   ! Batched through `trc_xc_grid_t`, which is the structure the XC
   ! integrator already uses for exactly this shape of sum: it hands back,
   ! per batch, the shells that reach it and a local function numbering, so
   ! the inner matrix is the batch's own few functions rather than the
   ! molecule's. Without that this is nao^2 per point and hopeless past a
   ! few atoms.
   !
   subroutine assemble(b, g, wv, h, error)
      type(trc_basis_t), intent(in) :: b
      type(dft_grid_t), intent(in) :: g
      real(dp), intent(in) :: wv(:)
      real(dp), intent(inout) :: h(b%nao, b%nao)
      type(error_t), intent(inout) :: error

      type(trc_xc_grid_t) :: xg
      real(dp), allocatable :: phi(:, :), phiw(:, :), hloc(:, :), wb(:)
      real(dp) :: chi(NCART_MAX), gchi(3, NCART_MAX), d(3), acc
      integer :: ib, p0, p1, npb, s0, s1, nloc, is, ish, np, nc, i, j, ip, l, mu, nu

      ! The batcher wants the weight it will carry; this one already has the
      ! potential folded in, and a point whose product is zero is no use to
      ! anyone, so it is passed as the weight and the screening it does on
      ! magnitude is screening on the right thing.
      call xg%build(g%n_points, g%coords, wv, b, SAP_BATCH, 1.0e-12_dp, 0.0_dp)

      do ib = 1, xg%nbatch
         p0 = xg%b_off(ib); p1 = xg%b_off(ib + 1) - 1
         npb = p1 - p0 + 1
         s0 = xg%b_shoff(ib); s1 = xg%b_shoff(ib + 1) - 1
         nloc = xg%b_aooff(ib + 1) - xg%b_aooff(ib)
         if (npb <= 0 .or. nloc <= 0) cycle
         allocate (phi(npb, nloc), phiw(npb, nloc), hloc(nloc, nloc), wb(npb))
         phi = 0.0_dp
         do ip = 1, npb
            wb(ip) = xg%w(p0 + ip - 1)
         end do
         do is = s0, s1
            ish = xg%b_sh(is)
            np = b%sh_np(ish)
            l = b%sh_l(ish)
            nc = ncart(l)
            do ip = 1, npb
               d = xg%r(:, p0 + ip - 1) - b%sh_r(:, ish)
               call shell_collocate(l, np, b%sh_e(1:np, ish), b%sh_c(1:np, ish), d, chi, gchi)
               do i = 1, nc
                  phi(ip, xg%b_shao(is) + i) = chi(i)
               end do
            end do
         end do
         !
         ! hloc = phi^T diag(wb) phi, THROUGH A GEMM.
         !
         ! Written out as a triple loop first, and that was the whole cost of
         ! the guess: nloc runs into the hundreds on a dense slab, so this is
         ! npts * nloc^2 and on the 123-atom silica slice it took the SAP
         ! guess to 188 s against SAD's 7.3. The arithmetic is a GEMM and
         ! BLAS does it two orders of magnitude faster; scaling one copy of
         ! phi by the weight keeps it a single call and sidesteps the sign
         ! split a syrk would need, since w*V is not positive.
         !
         do j = 1, nloc
            do ip = 1, npb
               phiw(ip, j) = wb(ip)*phi(ip, j)
            end do
         end do
         hloc = 0.0_dp
         call pic_gemm(phiw, phi, hloc, transa='T', transb='N')
         do j = 1, nloc
            nu = xg%b_ao(xg%b_aooff(ib) + j - 1)
            do i = 1, nloc
               mu = xg%b_ao(xg%b_aooff(ib) + i - 1)
               h(mu, nu) = h(mu, nu) + hloc(i, j)
            end do
         end do
         deallocate (phi, phiw, hloc, wb)
      end do
      call xg%release()
   end subroutine assemble

end module trc_sap
