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
! Water, 6-31G, RHF and PBE. The RHF energy is the one link_gfortran asserts
! (-75.983974469890), and PBE at level 3 the one the single-GPU driver gives
! (-76.2981056519); both must come out the same on any rank count, which is
! what makes a run under `mpiexec -np 4` a test and not a demonstration.
!
program scf_mpi
   use trc_boys, only: dp
   use pic_mpi_lib, only: comm_t, comm_world, pic_mpi_init, pic_mpi_finalize, allreduce, MPI_MAX, MPI_MIN
   use trc_api, only: trc_basis_t, trc_bind_device
   use trc_scf_driver, only: trc_scf_options_t, trc_scf_result_t, trc_scf_run
   use trc_test_basis, only: read_xyz, build_631g
   implicit none

   real(dp), parameter :: E_RHF = -75.983974469890_dp, E_PBE = -76.2981056519_dp
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
   integer(kind=8) :: cc, rate

   call pic_mpi_init()
   comm = comm_world()
   dev = trc_bind_device(comm%rank())

   ! Any xyz for a timing; only water carries reference energies.
   xyz = 'water.xyz'
   if (command_argument_count() >= 1) call get_command_argument(1, xyz)
   ! A second argument of any kind prints the iterations (rank 0 only).
   call read_xyz(trim(xyz), natm, zint, at_r)
   allocate (at_z(natm)); at_z = real(zint, dp)
   call build_631g(natm, zint, at_r, nsh, sh_l, sh_np, sh_e, sh_c, sh_r, maxnp)
   call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, natm, at_z, at_r, maxnp)
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
   call bas%release()
   if (comm%rank() == 0) print '(a)', merge("scf_mpi: PASS", "scf_mpi: FAIL", ok)
   call pic_mpi_finalize()
   if (.not. ok) stop 1
end program scf_mpi
