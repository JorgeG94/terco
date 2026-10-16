!
! The SCF over MPI ranks: the same energy on one rank, two and four.
!
! Every rank binds its own device, builds the same basis and grid, and runs
! the same iteration; the Fock build takes every nranks-th quartet item of
! the sorted work list from its rank and the result is summed across ranks
! after each build. So the density is identical on every rank, and the
! energies must be too -- bitwise, in principle, since the same numbers are
! summed in the same order; the check allows rounding.
!
! Water, 6-31G, RHF and PBE; both must come out the same on any rank count,
! which is what makes a run under `mpiexec -np 4` a test and not a
! demonstration.
!
! The basis loads through `trc_basis_from_json`, the one construction path
! -- E_RHF and E_PBE below are its own full-double-precision values (2026-
! 10-08), not `build_631g`'s old 7-significant-figure literal ones, which
! differed by 4.2e-9 Eh (RHF) and 8.6e-8 Eh (PBE, the DFT quadrature being
! the more sensitive of the two to it). The tolerance stays 1e-8: it is a
! correctness check, not a place to absorb a basis change by loosening it.
!
! A second case afterward: water/6-31G*, C-PCM (HF) and IEF-PCM (B3LYP),
! eps 78.3553. trc_pcm's pcm_operator splits its own way -- every rank
! keeps the whole cavity and K replicated (cheap next to the integrals)
! and only the two solute-surface integral passes are split by point
! range -- so this checks both the total energy and E_diel (res%e_pcm)
! against this program's own np 1 water/6-31G* values (2026-10-15).
!
program scf_mpi
   use pic_types, only: dp, int64
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init, pic_mpi_finalize, allreduce, MPI_MAX, MPI_MIN
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz, unpack_basis
   use trc_basis_json, only: trc_basis_from_json
   use trc_error, only: error_t
   implicit none

   real(dp), parameter :: E_RHF = -75.983974465718_dp
   !! trc_basis_from_json's own 6-31G water values, 2026-10-08.
   real(dp), parameter :: E_PBE = -76.298105565582_dp
   !! Same date, same basis path.
   real(dp), parameter :: PCM_EPS = 78.3553_dp
   real(dp), parameter :: E_HF_CPCM = -76.022097922513_dp
   !! water/6-31G*, this program's own np 1 values, 2026-10-15.
   real(dp), parameter :: ED_HF_CPCM = -0.012745024061_dp
   real(dp), parameter :: E_B3LYP_IEFPCM = -76.418997382217_dp
   real(dp), parameter :: ED_B3LYP_IEFPCM = -0.011333511625_dp
   type(comm_t) :: comm
   integer :: natm, nsh, maxnp, dev, ic
   integer, allocatable :: zint(:), sh_l(:), sh_np(:)
   real(dp), allocatable :: at_r(:, :), at_z(:), sh_e(:, :), sh_c(:, :), sh_r(:, :)
   type(trc_basis_t) :: bas
   type(trc_scf_options_t) :: opts
   type(trc_scf_result_t) :: res
   real(dp) :: e(2), emax, emin, ref
   character(len=8) :: names(2) = [character(len=8) :: "", "pbe"]
   logical :: ok
   character(len=256) :: xyz, xyz_arg
   real(dp) :: t0, t1
   integer(int64) :: cc, rate

   call pic_mpi_init()
   comm = comm_world()
   dev = trc_bind_device(comm%rank())

   ! Any xyz for a timing; only water carries reference energies.
   xyz = 'water.xyz'
   if (command_argument_count() >= 1) call get_command_argument(1, xyz)
   ! A second argument of any kind prints the iterations (rank 0 only).
   call read_xyz(trim(xyz), natm, zint, at_r)
   allocate (at_z(natm)); at_z = real(zint, dp)
   block
      type(error_t) :: err
      call trc_basis_from_json('basis_sets/6-31g.json', natm, zint, at_r, bas, err)
      if (err%has_error()) then
         print '(a)', 'scf_mpi: '//err%get_message(); stop 1
      end if
   end block
   call unpack_basis(bas, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   call bas%to_device()

   ok = .true.
   do ic = 1, 2
      opts = trc_scf_options_t()
      opts%functional = names(ic)
      opts%verbose = command_argument_count() >= 2
      ! Optional third and fourth arguments: level shift and damping weight.
      if (command_argument_count() >= 3) then
         call get_command_argument(3, xyz_arg); read (xyz_arg, *) opts%level_shift
      end if
      if (command_argument_count() >= 4) then
         call get_command_argument(4, xyz_arg); read (xyz_arg, *) opts%damp
      end if
      call system_clock(cc, rate); t0 = real(cc, dp)/real(rate, dp)
      call trc_scf_run(bas, nint(sum(at_z))/2, nint(sum(at_z))/2, opts, res, comm=comm)
      call system_clock(cc, rate); t1 = real(cc, dp)/real(rate, dp)
      e(ic) = res%energy
      emax = e(ic); emin = e(ic)
      ! pic-mpi's single-rank backend refuses a collective outright.
      if (comm%size() > 1) then
         call allreduce(comm, emax, op=MPI_MAX)
         call allreduce(comm, emin, op=MPI_MIN)
      end if
      ref = merge(E_RHF, E_PBE, ic == 1)
      if (comm%rank() == 0) then
         print '(a,i0,a,a5,a,f20.12,a,i0,a,es9.2,a,es9.2)', "scf_mpi: ", comm%size(), " rank(s)  ", &
            merge("rhf  ", names(ic)(1:5), ic == 1), "  E = ", e(ic), "  iter ", res%iterations, &
            "  spread over ranks ", emax - emin, "  vs reference ", abs(e(ic) - ref)
         print '(a,f8.2,a)', "         wall ", t1 - t0, " s"
      end if
      if (.not. res%converged .or. emax - emin > 1.0e-12_dp) ok = .false.
      if (trim(xyz) == 'water.xyz' .and. abs(e(ic) - ref) > 1.0e-8_dp) ok = .false.
   end do

   ! --- PCM: water/6-31G*, C-PCM (HF) and IEF-PCM (B3LYP) -----------------
   !
   ! pcm%operator (src/solvation/trc_pcm.F90) splits the two solute-surface integral
   ! passes over a contiguous range of points per rank; the cavity and K
   ! stay whole and replicated. Both the total energy AND E_diel (res%e_pcm)
   ! must land on this program's own np 1 values on every rank count.
   block
      type(trc_basis_t) :: bas_pcm
      type(trc_scf_options_t) :: opts_pcm
      type(trc_scf_result_t) :: res_pcm
      integer :: natm_pcm, jc
      integer, allocatable :: zint_pcm(:)
      real(dp), allocatable :: at_r_pcm(:, :), at_z_pcm(:)
      real(dp) :: ep(2), edp(2), emaxp, eminp, refe, refd
      character(len=8), parameter :: pcm_model(2) = [character(len=8) :: "cpcm", "iefpcm"]
      character(len=8), parameter :: pcm_fn(2) = [character(len=8) :: "", "b3lyp"]
      character(len=5), parameter :: pcm_tag(2) = [character(len=5) :: "cpcm ", "iefpc"]
      type(error_t) :: err_pcm

      call read_xyz('water.xyz', natm_pcm, zint_pcm, at_r_pcm)
      allocate (at_z_pcm(natm_pcm)); at_z_pcm = real(zint_pcm, dp)
      call trc_basis_from_json('basis_sets/6-31gs.json', natm_pcm, zint_pcm, at_r_pcm, bas_pcm, err_pcm)
      if (err_pcm%has_error()) then
         print '(a)', 'scf_mpi: '//err_pcm%get_message(); stop 1
      end if
      call bas_pcm%to_device()

      do jc = 1, 2
         opts_pcm = trc_scf_options_t()
         opts_pcm%functional = pcm_fn(jc)
         opts_pcm%pcm_model = pcm_model(jc)
         opts_pcm%pcm_eps = PCM_EPS
         call trc_scf_run(bas_pcm, nint(sum(at_z_pcm))/2, nint(sum(at_z_pcm))/2, opts_pcm, res_pcm, comm=comm)
         ep(jc) = res_pcm%energy
         edp(jc) = res_pcm%e_pcm
         emaxp = ep(jc); eminp = ep(jc)
         if (comm%size() > 1) then
            call allreduce(comm, emaxp, op=MPI_MAX)
            call allreduce(comm, eminp, op=MPI_MIN)
         end if
         refe = merge(E_HF_CPCM, E_B3LYP_IEFPCM, jc == 1)
         refd = merge(ED_HF_CPCM, ED_B3LYP_IEFPCM, jc == 1)
         if (comm%rank() == 0) then
            print '(a,i0,a,a5,a,f20.12,a,f16.12,a,i0,a,es9.2,a,es9.2)', "scf_mpi: ", comm%size(), &
               " rank(s)  ", pcm_tag(jc), "  E = ", ep(jc), "  E_diel = ", edp(jc), "  iter ", &
               res_pcm%iterations, "  spread over ranks ", emaxp - eminp, "  vs reference ", &
               abs(ep(jc) - refe)
         end if
         if (.not. res_pcm%converged .or. emaxp - eminp > 1.0e-12_dp) ok = .false.
         if (abs(ep(jc) - refe) > 1.0e-8_dp) ok = .false.
         if (abs(edp(jc) - refd) > 1.0e-8_dp) ok = .false.
      end do
      call bas_pcm%release()
   end block

   ! The same two runs through the C entry mqc uses, trc_scf_mpi, handed the
   ! world communicator's Fortran handle. It binds the device itself and
   ! must land on the energies above on every rank.
   block
      use, intrinsic :: iso_c_binding, only: c_int, c_double, c_ptr, c_null_ptr, c_null_char
      use trc_c_interfaces, only: trc_basis_create, trc_basis_destroy, trc_scf_mpi, TRC_OK
      integer(c_int) :: fcomm, rc, niter, nocc
      integer(c_int), allocatable :: lc(:), npc(:)
      real(c_double), allocatable :: zat(:), d_lib(:, :), eps_lib(:)
      real(c_double) :: e_lib, e_xc, e_c(2)
      type(c_ptr) :: hbas
      integer :: nao
      ! The handle comes from pic-mpi, not from `use mpi_f08`, so this
      ! builds on every backend -- including the legacy `mpi` module, which
      ! is all some MPIs (Cray MPICH under nvfortran) offer. On the serial
      ! backend it is a non-null placeholder, so the "collective" path is
      ! still exercised with no real MPI (rank 0 of 1).
      fcomm = int(comm%fortran_handle(), c_int)
      allocate (zat(natm), lc(nsh), npc(nsh))
      zat = real(at_z, c_double); lc = int(sh_l, c_int); npc = int(sh_np, c_int)
      rc = trc_basis_create(int(nsh, c_int), int(maxnp, c_int), lc, npc, sh_e, sh_c, sh_r, &
                            int(natm, c_int), zat, at_r, hbas)
      if (rc /= TRC_OK) ok = .false.
      nao = bas%nao
      nocc = int(nint(sum(at_z))/2, c_int)
      allocate (d_lib(nao, nao), eps_lib(nao))
      rc = trc_scf_mpi(fcomm, hbas, nocc, nocc, c_null_char, 3_c_int, 1.0e-10_c_double, 1.0e-7_c_double, &
                       100_c_int, c_null_ptr, 0_c_int, e_lib, e_xc, d_lib, eps_lib, niter)
      if (rc /= TRC_OK) ok = .false.
      e_c(1) = e_lib
      rc = trc_scf_mpi(fcomm, hbas, nocc, nocc, 'pbe'//c_null_char, 3_c_int, 1.0e-10_c_double, 1.0e-7_c_double, &
                       100_c_int, c_null_ptr, 0_c_int, e_lib, e_xc, d_lib, eps_lib, niter)
      if (rc /= TRC_OK) ok = .false.
      e_c(2) = e_lib
      rc = trc_basis_destroy(hbas)
      if (comm%rank() == 0) print '(a,es9.2,a,es9.2)', "scf_mpi: trc_scf_mpi vs driver  rhf ", &
         abs(e_c(1) - e(1)), "  pbe ", abs(e_c(2) - e(2))
      if (any(abs(e_c - e) > 1.0e-10_dp)) ok = .false.
      ! The same two runs through the context, split over the same ranks. GWH
      ! as the guess, which is what the driver above started from.
      block
         use trc_c_interfaces, only: trc_create, trc_destroy, trc_set_molecule, trc_set_basis_arrays, &
                                     trc_set_guess, trc_set_comm, trc_set_method, trc_run_scf, trc_energy, &
                                     TRC_GUESS_GWH
         type(c_ptr) :: h
         real(c_double) :: e_h(2)
         rc = trc_create(h)
         if (rc == TRC_OK) rc = trc_set_molecule(h, int(natm, c_int), zat, reshape(at_r, [3*natm]), 0_c_int, 1_c_int)
         if (rc == TRC_OK) rc = trc_set_basis_arrays(h, int(nsh, c_int), int(maxnp, c_int), lc, npc, &
                                                     reshape(sh_e, [size(sh_e)]), reshape(sh_c, [size(sh_c)]), &
                                                     reshape(sh_r, [size(sh_r)]))
         if (rc == TRC_OK) rc = trc_set_guess(h, TRC_GUESS_GWH, c_null_ptr, 1_c_int)
         if (rc == TRC_OK) rc = trc_set_comm(h, fcomm)
         if (rc == TRC_OK) rc = trc_run_scf(h)
         if (rc == TRC_OK) rc = trc_energy(h, e_h(1))
         if (rc == TRC_OK) rc = trc_set_method(h, 'pbe'//c_null_char, 3_c_int)
         if (rc == TRC_OK) rc = trc_run_scf(h)
         if (rc == TRC_OK) rc = trc_energy(h, e_h(2))
         if (rc /= TRC_OK) ok = .false.
         rc = trc_destroy(h)
         if (comm%rank() == 0) print '(a,es9.2,a,es9.2)', "scf_mpi: context vs driver      rhf ", &
            abs(e_h(1) - e(1)), "  pbe ", abs(e_h(2) - e(2))
         if (any(abs(e_h - e) > 1.0e-10_dp)) ok = .false.
      end block
   end block

   call bas%release()
   if (comm%rank() == 0) print '(a)', merge("scf_mpi: PASS", "scf_mpi: FAIL", ok)
   call pic_mpi_finalize()
   if (.not. ok) stop 1
end program scf_mpi
