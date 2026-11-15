module rccsd_test_support
   !!
   !! What the RI-CCSD(T) checks share (check_rccsd_ints, check_rccsd_inter): the
   !! two molecules and their reference files, the geometry reader that converts
   !! with CODATA 2018, the tight SCF, the MO dipole blocks and the sign-invariant
   !! probes of test/pyscf_refs/generate.py's gen_rccsd_ints.
   !!
   use pic_types, only: dp
   use trc_api, only: trc_basis_t, trc_pairlist_t, trc_multipoles, TRC_NMULT
   use trc_sph, only: sph_cart_to_sph
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_df_kernels, only: DF_LMAX, DF_LMAX_AUX
   use trc_constants, only: ANGSTROM_TO_BOHR
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   implicit none
   private
   public :: NMOL, MOL_LABEL, MOL_XYZ, MOL_REF, MOL_NFROZEN, TOL_PROBE, PROBE_FLOOR
   public :: read_xyz_bohr, load_molecule, run_scf, mo_dipoles, space_block, probe2, probe4

   integer, parameter :: NMOL = 2
   character(len=16), parameter :: MOL_LABEL(NMOL) = [character(len=16) :: "water", "h2o2 (C1)"]
   character(len=24), parameter :: MOL_XYZ(NMOL) = [character(len=24) :: "hirshfeld_water.xyz", "h2o2_c1.xyz"]
   character(len=24), parameter :: MOL_REF(NMOL) = [character(len=24) :: "rccsd_ints_w1.ref", "rccsd_ints_h2o2_c1.ref"]
   integer, parameter :: MOL_NFROZEN(NMOL) = [1, 2]

   real(dp), parameter :: TOL_PROBE = 1.0e-10_dp
   !! A probe matches when |terco - pyscf| <= TOL_PROBE max(|pyscf|, PROBE_FLOOR): relative, but a
   !! probe that is zero by symmetry is held to an absolute 1e-11 (the reference reproduces to 2e-12).
   real(dp), parameter :: PROBE_FLOOR = 0.1_dp

contains

   subroutine load_molecule(im, bas, aux, pl, zint, nocc, supported)
      !! Basis, RIFIT auxiliary basis and pair list of molecule `im`, on the device. `supported` is false
      !! (and nothing is built) when this build's DF kernels stop below the cc-pVDZ-RIFIT shells.
      integer, intent(in) :: im
      type(trc_basis_t), intent(out) :: bas, aux
      type(trc_pairlist_t), intent(out) :: pl
      integer, allocatable, intent(out) :: zint(:)
      integer, intent(out) :: nocc
      logical, intent(out) :: supported

      real(dp), allocatable :: at_r(:, :)
      integer :: natm
      type(error_t) :: err

      call read_xyz_bohr(trim(MOL_XYZ(im)), natm, zint, at_r)
      call trc_basis_from_json("basis_sets/cc-pvdz.json", natm, zint, at_r, bas, err)
      if (err%has_error()) error stop "rccsd_test_support: "//err%get_message()
      call trc_basis_from_json("basis_sets/cc-pvdz-rifit.json", natm, zint, at_r, aux, err)
      if (err%has_error()) error stop "rccsd_test_support: "//err%get_message()
      supported = maxval(bas%sh_l) <= DF_LMAX .and. maxval(aux%sh_l) <= DF_LMAX_AUX
      nocc = sum(zint)/2
      if (.not. supported) return
      call bas%to_device()
      call aux%to_device()
      call pl%build(bas, 1.0e-12_dp)
      call pl%to_device()
   end subroutine load_molecule

   subroutine run_scf(bas, nocc, res)
      !! RHF to energy 1e-13 and commutator 1e-11. The probes are first order in the orbital error (the
      !! energies second), so both sides converge the commutator to 1e-11: at 1e-9 the probes carried
      !! 1e-10 of SCF noise, at 1e-11 terco and pyscf agree to 3e-12.
      type(trc_basis_t), intent(in) :: bas
      integer, intent(in) :: nocc
      type(trc_scf_result_t), intent(out) :: res

      type(trc_scf_options_t) :: opts

      opts = trc_scf_options_t()
      opts%conv_energy = 1.0e-13_dp
      opts%conv_density = 1.0e-11_dp
      opts%conv_diis = 1.0e-11_dp
      call trc_scf_run(bas, nocc, nocc, opts, res)
   end subroutine run_scf

   function probe2(a, spaces, dip, no) result(p)
      !! |A|, sum A_pq X_pq, sum A_pq (XY)_pq (the reference probes of a 2-index tensor)
      real(dp), intent(in) :: a(:, :)
      character(len=2), intent(in) :: spaces
      real(dp), intent(in) :: dip(:, :, :)
      integer, intent(in) :: no
      real(dp) :: p(3)

      real(dp), allocatable :: xpq(:, :), xy(:, :)

      xpq = space_block(dip(:, :, 1), spaces(1:1), spaces(2:2), no)
      xy = space_block(dip(:, :, 4), spaces(1:1), spaces(2:2), no)
      p(1) = sqrt(sum(a*a))
      p(2) = sum(a*xpq)
      p(3) = sum(a*xy)
   end function probe2

   subroutine read_xyz_bohr(path, natm, zint, at_r)
      !! H and O only. Not trc_test_basis read_xyz: that helper converts with an older Bohr radius
      !! (1.8897261254578281 per Angstrom against 1.8897261246257702 for CODATA 2018, which the
      !! library and the pyscf reference use), a 4e-10 relative shift of the geometry that moves an
      !! SCF energy by 1e-10.
      character(len=*), intent(in) :: path
      integer, intent(out) :: natm
      integer, allocatable, intent(out) :: zint(:)
      real(dp), allocatable, intent(out) :: at_r(:, :)
      integer :: unit, i
      character(len=8) :: symbol
      character(len=256) :: line
      real(dp) :: xyz(3)

      open (newunit=unit, file=path, status="old", action="read")
      read (unit, *) natm
      read (unit, "(a)") line
      allocate (zint(natm), at_r(3, natm))
      do i = 1, natm
         read (unit, *) symbol, xyz
         select case (trim(symbol))
         case ("H"); zint(i) = 1
         case ("O"); zint(i) = 8
         case default; error stop "rccsd_test_support: only H and O are supported"
         end select
         at_r(:, i) = xyz*ANGSTROM_TO_BOHR
      end do
      close (unit)
   end subroutine read_xyz_bohr

   subroutine mo_dipoles(bas, pl, cact, dip)
      !! dip(:, :, 1:3) = C^T r C over the active orbitals (origin 0); dip(:, :, 4) = X Y
      type(trc_basis_t), intent(in) :: bas
      type(trc_pairlist_t), intent(in) :: pl
      real(dp), intent(in) :: cact(:, :)
      real(dp), intent(out) :: dip(:, :, :)

      real(dp), allocatable :: mm(:, :, :), msph(:, :, :), tmp(:, :)
      integer :: k, nact

      nact = size(cact, 2)
      allocate (mm(bas%nao, bas%nao, TRC_NMULT), msph(bas%nao_sph, bas%nao_sph, 3))
      !$acc enter data create(mm, msph)
      call trc_multipoles(bas, pl, [0.0_dp, 0.0_dp, 0.0_dp], mm)
      do k = 1, 3
         call sph_cart_to_sph(bas, mm(:, :, k), msph(:, :, k))
      end do
      !$acc update self(msph)
      !$acc exit data delete(mm, msph)
      allocate (tmp(bas%nao_sph, nact))
      do k = 1, 3
         tmp(:, :) = matmul(msph(:, :, k), cact)
         dip(:, :, k) = matmul(transpose(cact), tmp)
      end do
      dip(:, :, 4) = matmul(dip(:, :, 1), dip(:, :, 2))
   end subroutine mo_dipoles

   function space_block(m, sa, sb, no) result(blk)
      !! m restricted to the rows of space `sa` and the columns of `sb` (o = 1:no, v = no+1:)
      real(dp), intent(in) :: m(:, :)
      character, intent(in) :: sa, sb
      integer, intent(in) :: no
      real(dp), allocatable :: blk(:, :)
      integer :: r0, r1, c0, c1
      r0 = merge(1, no + 1, sa == "o")
      r1 = merge(no, size(m, 1), sa == "o")
      c0 = merge(1, no + 1, sb == "o")
      c1 = merge(no, size(m, 2), sb == "o")
      allocate (blk(r1 - r0 + 1, c1 - c0 + 1))
      blk = m(r0:r1, c0:c1)
   end function space_block

   function probe4(a, spaces, dip, no) result(p)
      !! |A|, sum A X_pq Y_rs, sum A X_pr Z_qs, sum A Y_ps Z_qr (the reference probes)
      real(dp), intent(in) :: a(:, :, :, :)
      character(len=4), intent(in) :: spaces
      real(dp), intent(in) :: dip(:, :, :)
      integer, intent(in) :: no
      real(dp) :: p(4)

      real(dp), allocatable :: xpq(:, :), yrs(:, :), xpr(:, :), zqs(:, :), yps(:, :), zqr(:, :)
      integer :: i1, i2, i3, i4
      real(dp) :: v

      xpq = space_block(dip(:, :, 1), spaces(1:1), spaces(2:2), no)
      yrs = space_block(dip(:, :, 2), spaces(3:3), spaces(4:4), no)
      xpr = space_block(dip(:, :, 1), spaces(1:1), spaces(3:3), no)
      zqs = space_block(dip(:, :, 3), spaces(2:2), spaces(4:4), no)
      yps = space_block(dip(:, :, 2), spaces(1:1), spaces(4:4), no)
      zqr = space_block(dip(:, :, 3), spaces(2:2), spaces(3:3), no)
      p = 0.0_dp
      do i4 = 1, size(a, 4)
         do i3 = 1, size(a, 3)
            do i2 = 1, size(a, 2)
               do i1 = 1, size(a, 1)
                  v = a(i1, i2, i3, i4)
                  p(1) = p(1) + v*v
                  p(2) = p(2) + v*xpq(i1, i2)*yrs(i3, i4)
                  p(3) = p(3) + v*xpr(i1, i3)*zqs(i2, i4)
                  p(4) = p(4) + v*yps(i1, i4)*zqr(i2, i3)
               end do
            end do
         end do
      end do
      p(1) = sqrt(p(1))
   end function probe4

end module rccsd_test_support
