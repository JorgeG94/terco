program trc_main
   !! The `terco` command: an SCF on an xyz file, from the shell.
   !!
   !! `terco molecule.xyz [options]`, in the spirit of xtb's command line.
   !!
   !! Everything below the argument handling is the library's context API
   !! (trc_capi): the guesses, the SCF, RI-MP2 and the energies are exactly
   !! what a C or Python host gets. `terco --help` lists the options.
   !!
   !! This program itself has NO internal (`contains`) procedures -- see
   !! trc_run's header comment for why that is load-bearing, not a style
   !! choice, the moment more than one rank is involved.
   use pic_mpi_lib, only: comm_world, pic_mpi_init
   use pic_profiler, only: profiler_init
   use trc_api, only: trc_bind_device_env
   use trc_run, only: run_t, run_define_cli, run_parse, run_calculation, run_quit
   implicit none

   type(run_t) :: r
   integer :: dev

   ! MPI, in every build: pic-mpi's single-rank backend (TERCO_ENABLE_MPI=OFF)
   ! makes this a no-op, so the rest of the program need not branch on it.
   ! `r%wcomm` and `r%nranks` answer everything the rest of the program needs
   ! to know about ranks; `r%rank0` is what gates every print so N
   ! ranks produce ONE copy of the output, not N.
   !
   ! The device first, from the launcher's node-local rank: a GPU-aware MPI
   ! records the current device in MPI_Init, and switching afterwards costs
   ! Cray MPICH its CUDA IPC (trc_bind_device_env). The split-communicator
   ! bind in run_calculate_setup still runs, and is a no-op when it lands on this device.
   dev = trc_bind_device_env()
   call pic_mpi_init()
   ! pic's region profiler, which is what puts the library's named NVTX
   ! ranges on an Nsight Systems timeline (a TERCO_NVTX build); without NVTX
   ! it only keeps host wall times nobody reads, so it is always on.
   call profiler_init()
   r%wcomm = comm_world()
   r%nranks = r%wcomm%size()
   r%rank0 = r%wcomm%leader()

   call r%wall%start()

   call run_define_cli(r)
   call run_parse(r)
   call run_calculation(r)
   call run_quit(r, 0)

end program trc_main
