module trc_sap
   !! The superposition-of-atomic-potentials guess
   !!
   !! SAP builds a one-electron Hamiltonian from free-atom potentials and
   !! diagonalises it, rather than superposing densities and diagonalising the
   !! Fock matrix they make. Lehtola (JCTC 15 1593, 2019): the atomic potential is
   !! what an electron in the molecule sees, and carries the atomic shell
   !! structure into the guess orbitals directly.
   !!
   !!   H_sap = T + V_ne + sum_A Vscr_A(|r - R_A|)
   !!
   !! Two qualifications. First, this is not the cheap SAP: Lehtola's uses
   !! tabulated potentials with no atomic calculation. No table is in this tree
   !! (and inventing one is not to be done quietly), so `Vscr_A` comes from the
   !! free-atom SCF SAD already runs and caches, and SAP costs SAD plus a
   !! quadrature: it is for guess quality, not speed. A real table would change
   !! only `atom_screen`.
   !!
   !! Second, a Hartree-Fock atom's effective potential is nonlocal; this uses a
   !! local approximation (the atom's Hartree potential plus the Dirac-Slater
   !! exchange of its density), as the tabulated potentials do. That is why SAP is
   !! a guess and not a method.
   !!
   !! The nuclear attraction -Z/r is singular and already exact, so it is not
   !! quadratured: T and V_ne come in as analytic matrices and only the smooth,
   !! short-ranged screening potential is integrated.
   use pic_types, only: dp, int64
   use pic_logger, only: global_logger
   use pic_timer, only: timer_type
   use trc_constants, only: PI
   use trc_api, only: trc_basis_t
   use trc_error, only: error_t, ERROR_VALIDATION
   use trc_sad, only: trc_sad_build, trc_sadq_build, trc_atom_ranges
   use trc_dft_radial, only: dft_radial_treutler_ahlrichs
   use trc_lebedev, only: lebedev_grid, lebedev_order_at_least
   use trc_dft_grid, only: dft_grid_t, dft_grid_build
   use trc_xc_batch, only: trc_xc_grid_t
   use trc_collocation, only: collocation_shell, NCART_MAX
   use pic_blas_interfaces, only: pic_gemm
   use trc_api, only: trc_pairlist_t, trc_df_3c
   use trc_basis_json, only: trc_fitted_shell_from_json
   implicit none
   private

   public :: trc_sap_build, trc_sap_grid, trc_sap_screen_table, SAP_BASIS_DEFAULT

   integer, parameter :: SAP_NRAD = 400
   !! Radial nodes per element for the screening table: generous rather than
   !! tuned (the potential is smooth, M4 puts nodes where it bends); built once.

   integer, parameter :: SAP_NANG = 26
   !! Angular order for the spherical average of the atomic density. The
   !! free-atom density is already spherical; averaging is free and makes a
   !! non-spherical atom degrade the guess rather than silently tilt it.

   character(len=*), parameter :: SAP_BASIS_DEFAULT = "basis_sets/sap_helfem_large.json"
   !! The fitted potential set, unless a caller names another: Lehtola's fully
   !! numerical HelFEM fit (the GRASP sets are relativistic, no gain below heavy elements).

   integer(int64), parameter :: SAP_TENSOR_BYTES = 2147483648_int64
   !! Device bytes the three-centre tensor may take at once.
   !!
   !! `trc_df_3c` returns the whole (nao, nao, naux) tensor with naux one function
   !! per atom, so all atoms at once wants nao^2 * natm doubles (17 GB for 150
   !! waters in cc-pVDZ). The sum over atoms is order-free, so it is grouped to
   !! fit here.
   !!
   !! A memory budget only: pass count does not matter (silica slice, sixteen
   !! passes to one: 21.78 s against 21.77), since the three-centre integrals are
   !! the whole cost. Set well under a card because the SCF's resident arrays
   !! are already on it.

   integer, parameter :: SAP_BATCH = 128
   !! Points per quadrature batch for the assembly.

   !
   integer, parameter :: SAP_GRID_LEVEL = 3
   !! Quadrature level for the screening integral.
   !!
   !! This is the DFT default, deliberately. A cheaper grid was tried first, on
   !! the reasoning that a throwaway guess needs no converged XC grid, but
   !! quadrature error here is a wrong operator, not a small perturbation. On the
   !! 123-atom silica slice, cc-pVDZ, at iteration 12:
   !!
   !!              E(12)                 dE(12)     guess cost
   !!   level 0    -11194.245761258     -9.5e-05      7 s
   !!   level 1    -11194.245769153     -9.9e-05     25 s
   !!   level 3    -11194.245774662     -3.3e-10     80 s
   !!   SAD        -11194.245774661     -6.0e-09      <1 s
   !!
   !! Level 3 beats SAD by about one iteration; levels 0 and 1 are not converging.
   !! So SAP is currently a loss at this size (it saves ~12 s and costs 80, of
   !! which 18 is `assemble` and 4 the per-point potential, both host-side and
   !! single-threaded). The fix is to move those to the device as trc_xc does, not
   !! to coarsen the grid. On small molecules it earns its place (water/6-31G: 9
   !! iterations against 10, no measurable cost).
   !!
   !! TRC_SAP_LEVEL overrides it, and an explicit `level` argument overrides both.

contains

   pure integer function ncart(l) result(res)
      integer, intent(in) :: l
      res = (l + 1)*(l + 2)/2
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
   ! `nel` comes back as 4 pi int r^2 rho dr: a check on the quadrature and the
   ! coefficient of the tail, since past the last node the atom is a point charge
   ! -nel, so Vscr = nel/r and V_ne + Vscr goes to zero for a neutral atom.
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
      real(dp) :: s, acc_in, acc_out, rh
      real(dp) :: pt(3), d(3)
      integer :: nrad, nang, k, ia, ish, np, nc, off, nz, i, j

      nz = size(dblk, 1)
      call dft_radial_treutler_ahlrichs(SAP_NRAD, z, r, dr, error)
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
               call collocation_shell(b%sh_l(ish), np, b%sh_e(1:np, ish), b%sh_c(1:np, ish), d, chi, gchi)
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
   ! Vscr of one element at an arbitrary distance: linear between nodes, nel/rr
   ! past the last one (a point charge), and the first node's value inside it
   ! (unresolved by any molecular grid; the nuclear term dominates there).
   !
   pure real(dp) function interp(r, v, nel, rr) result(res)
      real(dp), intent(in) :: nel, rr
      real(dp), intent(in) :: r(:), v(:)
      integer :: lo, hi, mid
      real(dp) :: t
      if (rr <= r(1)) then
         res = v(1)
      else if (rr >= r(size(r))) then
         res = nel/rr
      else
         lo = 1
         hi = size(r)
         do while (hi - lo > 1)
            mid = (lo + hi)/2
            if (r(mid) > rr) then
               hi = mid
            else
               lo = mid
            end if
         end do
         t = (rr - r(lo))/(r(hi) - r(lo))
         res = v(lo)*(1.0_dp - t) + v(hi)*t
      end if
   end function interp

   !
   ! The screening potential from FITTED atomic charge distributions.
   !
   ! Lehtola's Gaussian fits make this a three-centre contraction instead of a
   ! quadrature:
   !
   !   Vscr(mu,nu) = sum_A (mu nu | rho_A),   rho_A = sum_p c_p g_p
   !
   ! with one contracted s shell per atom (naux = atom count). It is the same
   ! integral density fitting does, on the device, and replaces a host-side grid
   ! assembly that cost 80 s on a 123-atom slab; the grid path below is kept as
   ! the independent check.
   !
   ! NORMALISATION, THE TRAP. The published coefficients are charges over
   ! unit-normalised distributions, summing to -Z (checked by
   ! trc_fitted_shell_from_json). Two things must happen:
   !
   ! Sign: `(mu nu | P)` is for a POSITIVE unit charge, but an electron in the
   ! field of Z electrons has energy +Z/R, so the coefficient passed is
   ! -c_published and the far field is +Z/R.
   !
   ! Normalisation: `basis_build` multiplies by common_fac_sp(l) only, so the
   ! engine integrates (sh_c / common_fac_sp) exp(-alpha r^2), and a unit charge
   ! distribution is (alpha/pi)^(3/2) exp(-alpha r^2). The coefficient passed is
   ! therefore -c_published * (alpha/pi)^(3/2), NOT the wavefunction
   ! normalisation (2 alpha/pi)^(3/4) used for every other basis, which gives a
   ! smooth exponent-dependent error that still converges, silently.
   !
   subroutine trc_sap_build(b, pl, tmat, vmat, hsap, error, verbose, path)
      type(trc_basis_t),    intent(in) :: b
      type(trc_pairlist_t), intent(in) :: pl
      real(dp), intent(in)  :: tmat(b%nao, b%nao), vmat(b%nao, b%nao)
      !! The analytic one-electron matrices; only the screening potential is
      !! computed here.
      real(dp), intent(out) :: hsap(b%nao, b%nao)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose
      character(len=*), intent(in), optional :: path

      type(trc_basis_t) :: aux
      real(dp), allocatable :: e(:), c(:), tens(:, :, :)
      real(dp), allocatable :: a_e(:, :), a_c(:, :), a_r(:, :)
      integer,  allocatable :: a_l(:), a_np(:)
      character(len=:), allocatable :: bpath
      integer :: ia, i, j, k, np, nmax, n0, n1, nch, ish, nao, naux, chunk
      real(dp) :: cfac, acc
      logical :: talk

      talk = .false.
      if (present(verbose)) talk = verbose
      nao = b%nao
      bpath = SAP_BASIS_DEFAULT
      block
         character(len=512) :: pe
         pe = " "
         call get_environment_variable("TRC_SAP_BASIS", pe)
         if (len_trim(pe) > 0) bpath = trim(pe)
      end block
      if (present(path)) bpath = path
      hsap = 0.0_dp
      !$acc enter data copyin(hsap)

      ! Widest fit in the molecule, so one rectangular block holds any chunk.
      ! A ghost (b%at_z(ia) == 0) has no nucleus and so no fitted potential
      ! to look up -- it is given a zero-coefficient stand-in below instead.
      nmax = 0
      do ia = 1, b%natm
         if (b%at_z(ia) == 0.0_dp) cycle
         call trc_fitted_shell_from_json(bpath, b%at_elem(ia), e, c, error)
         if (error%has_error()) return
         nmax = max(nmax, size(e))
      end do

      !
      ! What the builder does to a coefficient, asked rather than assumed:
      ! `basis_build` folds common_fac_sp(l) in, so one probe shell with
      ! coefficient 1 is built and the factor read back, instead of restating a
      ! constant from trc_api that could go stale (as trc_decontract does too).
      ! Worth it because the error is invisible: dividing instead of multiplying
      ! made the potential a uniform 2*sqrt(pi) too weak and the SCF still
      ! converged to the same energy; only the grid construction caught it.
      !
      block
         real(dp) :: p_e(1, 1), p_c(1, 1), p_r(3, 1), p_z(1)
         integer :: p_l(1), p_np(1)
         type(trc_basis_t) :: probe
         p_l = 0
         p_np = 1
         p_e = 1.0_dp
         p_c = 1.0_dp
         p_r = 0.0_dp
         p_z = 1.0_dp
         call probe%build(1, p_l, p_np, p_e, p_c, p_r, 1, p_z, p_r, 1)
         cfac = probe%sh_c(1, 1)
         call probe%release()
      end block
      if (cfac == 0.0_dp) then
         call error%set(ERROR_VALIDATION, "trc_sap: the basis builder zeroed a unit coefficient")
         return
      end if
      if (talk) then
         block
            character(len=1024) :: logline
            write (logline, "(a,a,a,i0,a)") "  sap: fitted potentials from ", trim(bpath), &
               ", up to ", nmax, " primitives per atom"
            call global_logger%info(trim(logline))
         end block
      end if

      chunk = int(max(1_int64, SAP_TENSOR_BYTES/(int(nao, int64)*int(nao, int64)*8_int64)))
      chunk = min(chunk, b%natm)
      if (talk) then
         block
            character(len=1024) :: logline
            write (logline, "(a,i0,a,i0,a)") "  sap: ", (b%natm + chunk - 1)/chunk, &
               " three-centre pass(es), ", chunk, " atoms each"
            call global_logger%info(trim(logline))
         end block
      end if
      do n0 = 1, b%natm, chunk
         n1 = min(n0 + chunk - 1, b%natm)
         nch = n1 - n0 + 1
         allocate (a_l(nch), a_np(nch), a_e(nmax, nch), a_c(nmax, nch), a_r(3, nch))
         a_l = 0
         a_e = 0.0_dp
         a_c = 0.0_dp
         do ish = 1, nch
            ia = n0 + ish - 1
            a_r(:, ish) = b%at_r(:, ia)
            if (b%at_z(ia) == 0.0_dp) then
               ! A ghost: no nucleus, so no screening potential -- a single
               ! zero-coefficient primitive contributes exactly nothing to
               ! `tens` below, which is what "SAP must skip ghost
               ! potentials" means here.
               a_np(ish) = 1
               a_e(1, ish) = 1.0_dp
               a_c(1, ish) = 0.0_dp
               cycle
            end if
            call trc_fitted_shell_from_json(bpath, b%at_elem(ia), e, c, error)
            if (error%has_error()) then
               deallocate (a_l, a_np, a_e, a_c, a_r)
               return
            end if
            np = size(e)
            a_np(ish) = np
            do k = 1, np
               a_e(k, ish) = e(k)
               a_c(k, ish) = -c(k)*(e(k)/PI)**1.5_dp/cfac
            end do
         end do
         ! The aux centres are the real ones; the builder ignores the nucleus.
         call aux%build(nch, a_l, a_np, a_e, a_c, a_r, nch, b%at_z(n0:n1), a_r, nmax)
         call aux%to_device()
         allocate (tens(b%nao, b%nao, aux%nao))
         !$acc enter data create(tens)
         call trc_df_3c(b, pl, aux, tens)
         ! Reduced where it was computed: only the sum over atoms is wanted, so
         ! only hsap crosses the bus, once, at the end.
         naux = aux%nao
         do concurrent(j=1:nao, i=1:nao) local(acc, k)
            acc = 0.0_dp
            do k = 1, naux
               acc = acc + tens(i, j, k)
            end do
            hsap(i, j) = hsap(i, j) + acc
         end do
         !$acc exit data delete(tens)
         deallocate (tens, a_l, a_np, a_e, a_c, a_r)
         call aux%release()
      end do
      !$acc update self(hsap)
      !$acc exit data delete(hsap)
      hsap = hsap + tmat + vmat
   end subroutine trc_sap_build

   !
   ! H_sap = T + V_ne + sum_A Vscr_A.
   !
   ! `tmat` and `vmat` are the analytic matrices from trc_1e. `nelec` present
   ! selects SADQ atoms (charge spread over the atoms) rather than neutral SAD ones.
   !
   subroutine trc_sap_grid(b, tmat, vmat, hsap, error, verbose, nelec, level)
      type(trc_basis_t), intent(in) :: b
      real(dp), intent(in) :: tmat(b%nao, b%nao), vmat(b%nao, b%nao)
      real(dp), intent(out) :: hsap(b%nao, b%nao)
      type(error_t), intent(inout) :: error
      logical, intent(in), optional :: verbose
      integer, intent(in), optional :: nelec   !! molecular electron count, for SADQ atoms
      integer, intent(in), optional :: level   !! quadrature level, as dft_grid_build spells it

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
      type(timer_type) :: tm
      logical :: talk

      talk = .false.
      if (present(verbose)) talk = verbose
      hsap = 0.0_dp
      call tm%start()

      ! --- the free atoms, and their ranges -------------------------------
      if (present(nelec)) then
         call trc_sadq_build(b, nelec, dguess, error, verbose=verbose)
      else
         call trc_sad_build(b, dguess, error, verbose=verbose)
      end if
      if (error%has_error()) return
      call trc_atom_ranges(b, atom_of, first, last, error)
      if (error%has_error()) return
      call stage("free atoms   ")

      ! --- one screening table per distinct (element, function count) ------
      allocate (tab(b%natm), tab_of(b%natm), zint(b%natm))
      tab_of = 0
      ntab = 0
      do ia = 1, b%natm
         if (first(ia) == 0) cycle
         if (b%at_z(ia) == 0.0_dp) cycle   ! a ghost: no nucleus, no screening potential
         zint(ia) = b%at_elem(ia)
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
            tab(ntab)%z = zint(ia)
            tab(ntab)%nz = nz
            if (talk) then
               block
                  character(len=1024) :: logline
                  write (logline, "(a,i0,a,f12.8,a,f12.8)") "  sap: Z = ", zint(ia), &
                     "  radial density integrates to ", tab(ntab)%nel, &
                     "   r*Vscr at the last node ", tab(ntab)%r(size(tab(ntab)%r))* &
                     tab(ntab)%v(size(tab(ntab)%v))
                  call global_logger%info(trim(logline))
               end block
            end if
            tab_of(ia) = ntab
         end if
      end do

      call stage("screen tables")

      ! --- the molecular quadrature ---------------------------------------
      lv = SAP_GRID_LEVEL
      if (present(level)) lv = level
      block
         character(len=16) :: e
         integer :: ios, lv_env
         e = " "
         call get_environment_variable("TRC_SAP_LEVEL", e)
         if (len_trim(e) > 0 .and. .not. present(level)) then
            read (e, *, iostat=ios) lv_env
            if (ios == 0) lv = lv_env
         end if
      end block
      ! Grid SIZE follows the element (ghosts keep their quadrature points);
      ! `zint` above, the real nuclear charge, is only the screening table's
      ! own key and stays 0/unset for a ghost.
      call dft_grid_build(b%at_r(:, 1:b%natm), b%at_elem(1:b%natm), g, error, level=lv)
      if (error%has_error()) return
      if (talk) then
         block
            character(len=1024) :: logline
            write (logline, "(a,i0,a,i0,a)") "  sap: quadrature level ", lv, ", ", &
               g%n_points, " points"
            call global_logger%info(trim(logline))
         end block
      end if

      ! Every atom contributes to every point and nothing is screened: Vscr falls
      ! off as nel/r, which is the long-range neutrality of V_ne + Vscr, and
      ! dropping it would leave a spurious charge.
      call stage("grid         ")
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

      call stage("potential    ")
      call assemble(b, g, vpt, hsap, error)
      if (error%has_error()) return
      call stage("assemble     ")

      hsap = hsap + tmat + vmat
      call g%destroy()
      deallocate (dguess, atom_of, first, last, zint, tab, tab_of, vpt)

   contains

      subroutine stage(what)
         !! Stage timing, printed only when TRC_BUILD_TIMING is set, as the Fock
         !! build's own stages are. Four routines with very different costs and
         !! no way to tell from the outside which one is slow.
         character(len=*), intent(in) :: what
         character(len=8) :: e
         e = " "
         call get_environment_variable("TRC_BUILD_TIMING", e)
         if (len_trim(e) > 0) then
            block
               character(len=1024) :: logline
               write (logline, "(a,a,f9.3,a)") "  [sap] ", what, &
                  tm%get_elapsed_time(), " s"
               call global_logger%info(trim(logline))
            end block
         end if
         call tm%start()
      end subroutine stage
   end subroutine trc_sap_grid

   !
   ! H_munu += sum_g wv(g) chi_mu(g) chi_nu(g).
   !
   ! Batched through `trc_xc_grid_t`, as the XC integrator does: per batch it
   ! returns the shells that reach it and a local numbering, so the inner matrix
   ! is the batch's few functions, not nao^2 per point.
   !
   subroutine assemble(b, g, wv, h, error)
      type(trc_basis_t), intent(in) :: b
      type(dft_grid_t), intent(in) :: g
      real(dp), intent(in) :: wv(:)
      real(dp), intent(inout) :: h(b%nao, b%nao)
      type(error_t), intent(inout) :: error

      type(trc_xc_grid_t) :: xg
      real(dp), allocatable :: phi(:, :), phiw(:, :), hloc(:, :), wb(:)
      real(dp) :: acc
      real(dp) :: chi(NCART_MAX), gchi(3, NCART_MAX), d(3)
      integer :: ib, p0, p1, npb, s0, s1, nloc, is, ish, np, nc, i, j, ip, l, mu, nu

      ! The batcher screens on the weight it carries, so pass the product with the
      ! potential already folded in.
      call xg%build(g%n_points, g%coords, wv, b, SAP_BATCH, 1.0e-12_dp, 0.0_dp)

      do ib = 1, xg%nbatch
         p0 = xg%b_off(ib)
         p1 = xg%b_off(ib + 1) - 1
         npb = p1 - p0 + 1
         s0 = xg%b_shoff(ib)
         s1 = xg%b_shoff(ib + 1) - 1
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
               call collocation_shell(l, np, b%sh_e(1:np, ish), b%sh_c(1:np, ish), d, chi, gchi)
               do i = 1, nc
                  phi(ip, xg%b_shao(is) + i) = chi(i)
               end do
            end do
         end do
         !
         ! hloc = phi^T diag(wb) phi, through a GEMM. A triple loop was the whole cost
         ! of the guess (nloc runs into the hundreds on a dense slab: 188 s against
         ! SAD's 7.3 on the 123-atom silica slice). Scaling one copy of phi by the
         ! weight keeps it a single call and avoids the sign split a syrk would
         ! need, since w*V is not positive.
         !
         do j = 1, nloc
            do ip = 1, npb
               phiw(ip, j) = wb(ip)*phi(ip, j)
            end do
         end do
         hloc = 0.0_dp
         call pic_gemm(phiw, phi, hloc, transa="T", transb="N")
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
