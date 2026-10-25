!
! Whether MPI is handed device buffers directly or through the host.
!
! A CUDA-aware MPI (NVHPC's HPC-X OpenMPI is one) takes a device address
! and moves the data itself, over NVLink between GPUs that have it; any
! other MPI has to be given host memory, so a device array goes through
! `!$acc update self`, MPI and `!$acc update device`. On two V100s joined by
! NVLink, RI-MP2's 2 GB block of B took 2.85 s through the host.
!
! The choice is made at run time so that one build runs on any MPI. The
! default is the build's (CMake option TERCO_MPI_DEVICE, on for GPU builds
! with MPI, which defines TRC_MPI_DEVICE), and the environment variable
! TERCO_MPI_DEVICE=0 or 1 overrides it. A host build has no device buffers
! and always answers no. An MPI that is not CUDA-aware, given a device
! address, crashes (Cray MPICH needs MPICH_GPU_SUPPORT_ENABLED=1 first), so
! the variable is the way out on such a machine.
!
module trc_mpi_device
   implicit none
   private

   public :: trc_mpi_device_direct

contains

   function trc_mpi_device_direct() result(direct)
      logical :: direct
      character(len=8) :: e
      integer :: n, st
#if defined(TRC_DEVICE) && TRC_MPI_DEVICE
      direct = .true.
#else
      direct = .false.
#endif
#ifdef TRC_DEVICE
      call get_environment_variable("TERCO_MPI_DEVICE", e, n, st)
      if (st == 0 .and. n > 0) direct = e(1:1) /= "0"
#else
      e = ""
      n = 0
      st = 0
#endif
   end function trc_mpi_device_direct

end module trc_mpi_device
