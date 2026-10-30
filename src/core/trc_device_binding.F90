!
! Which GPU this process uses: one device per rank, chosen from the rank
! within its node (trc_bind_device), or from the launcher's node-local rank
! variable before MPI is up (trc_bind_device_env). A host build has nothing
! to bind and every function returns -1.
!
! A core module so that anything may bind without reaching up into the API;
! `trc_api` re-exports all three functions.
!
module trc_device_binding
   implicit none
   private
   public :: trc_bind_device, trc_bind_device_env, trc_local_rank_env

   integer :: bound_dev = -1
   !! The device trc_bind_device last made current, -1 before any. A later
   !! call that picks the same device leaves it alone (see trc_bind_device_env).

contains

   !
   ! One device per rank: rank r takes device r modulo the number of
   ! devices. `rank` must be the caller's rank WITHIN ITS NODE, not
   ! necessarily its rank in the whole communicator: a communicator is not
   ! guaranteed to be node-ordered (a DDI compute communicator or a
   ! fragment's sub-communicator need not be), so a caller whose
   ! communicator might span more than one node should split it by shared
   ! memory first (pic-mpi's comm_t%split(), MPI_COMM_TYPE_SHARED
   ! underneath) and pass the rank in THAT communicator. trc_capi's
   ! trc_build_comm does exactly this for every fcomm-taking C entry.
   ! Call it right after MPI is up and BEFORE any basis, pair list
   ! or grid goes to the device -- those land on whichever device is current
   ! at the time. Under a host build there is nothing to bind. Returns the
   ! device chosen, -1 on a host build.
   !
   function trc_bind_device(rank) result(dev)
#ifdef TRC_DEVICE
      use openacc, only: acc_set_device_num, acc_get_num_devices, acc_device_nvidia
#endif
      integer, intent(in) :: rank
      integer :: dev
      dev = -1
#ifdef TRC_DEVICE
      block
         integer :: ndev
         ndev = acc_get_num_devices(acc_device_nvidia)
         if (ndev > 0) then
            dev = mod(rank, ndev)
            if (dev /= bound_dev) call acc_set_device_num(dev, acc_device_nvidia)
            bound_dev = dev
         end if
      end block
#endif
   end function trc_bind_device

   !
   ! trc_bind_device from the node-local rank the launcher exports, for use
   ! BEFORE MPI_Init. A GPU-aware MPI records the current device when it
   ! starts. Cray MPICH's GTL then warns "process is not using the same
   ! device it did during gtl_init" and disables CUDA IPC between the
   ! node's ranks when the process switches devices after MPI_Init, which
   ! is what binding from a split communicator does. Returns -1 and binds
   ! nothing when no launcher variable is set or on a host build. Once MPI
   ! is up the caller still binds from the split communicator as before.
   ! That call is a no-op when it picks the same device. If it picks a
   ! different one (a variable inherited from an outer launcher, say),
   ! it rebinds, so the result is never wrong, only slower to communicate.
   !
   function trc_bind_device_env() result(dev)
      integer :: dev, lr
      dev = -1
      lr = trc_local_rank_env()
      if (lr >= 0) dev = trc_bind_device(lr)
   end function trc_bind_device_env

   !
   ! This process's rank within its node, as the launcher exports it, -1
   ! when none does. The launcher-specific variables come first and
   ! SLURM_LOCALID last. An mpirun started inside a Slurm step inherits
   ! that step's SLURM_LOCALID (0 on every rank) but sets its own
   ! variable, so checking Slurm first would put every rank on device 0.
   !
   function trc_local_rank_env() result(lr)
      integer :: lr
      character(len=*), parameter :: VARS(*) = [character(len=26) :: &
                                                "OMPI_COMM_WORLD_LOCAL_RANK", "PALS_LOCAL_RANKID", &
                                                "MPI_LOCALRANKID", "MV2_COMM_WORLD_LOCAL_RANK", "SLURM_LOCALID"]
      character(len=32) :: value
      integer :: i, n, st, ios
      lr = -1
      do i = 1, size(VARS)
         call get_environment_variable(trim(VARS(i)), value, n, st)
         if (st /= 0 .or. n == 0) cycle
         read (value(1:n), *, iostat=ios) lr
         if (ios == 0 .and. lr >= 0) return
         lr = -1
      end do
   end function trc_local_rank_env

end module trc_device_binding
