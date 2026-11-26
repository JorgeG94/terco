module trc_device_memory
   !!
   !! Does an allocation fit on the device now? One place for the free-memory query and the refusal text, so a
   !! module that is about to map a large array can refuse (a status) where the runtime would abort.
   !!
   !! The free memory is the OpenACC runtime's `acc_property_free_memory`, read at the time of the call; on a host
   !! build, or when the query gives nothing, it is -1 and every request fits. A caller asks about its OWN bytes
   !! before it maps them, and about the largest set of arrays it holds at once if it maps and frees in turn.
   !!
   !! POOL. The nvfortran runtime keeps memory it has freed in a pool. No margin is added for that here: the figure
   !! is exactly what the query returns. Measured (26.5, one 32 GB V100): naphthalene cc-pVTZ maps its (T) scratch
   !! of about 20 GB right after the CCSD work arrays are released and is not refused, so freed memory is reported
   !! free; caffeine cc-pVTZ, which aborted in CUDA before these checks, is refused at the W scratch.
   !!
   !! TEST HOOK. `device_memory_set_cap(n)` with n > 0 makes the query answer n instead of asking the device, so a
   !! host build can exercise every refusal. The RI-CCSD options (`trc_rccsd_options_t%device_free_cap`, and the
   !! triples options' `device_free_cap` of `trc_rccsd_triples`) set it for the duration of a call; 0 (the default)
   !! restores the real query. Not for production use: the value is never read from the environment.
   !!
   !! LOCAL MEMORY. Before it launches a kernel, CUDA holds one stack frame for every thread that can be resident on
   !! the device (multiprocessors times resident threads each), sized by the largest frame launched so far, and
   !! keeps it until `cudaDeviceSetLimit` shrinks it. A launch that cannot get it aborts, so a caller that checks
   !! before a kernel with a large frame also asks `device_memory_local_growth` what is still to be taken.
   use, intrinsic :: iso_c_binding, only: c_int, c_size_t
   use pic_types, only: dp, int64
   use pic_format, only: to_string_fixed
   use trc_error, only: error_t, ERROR_VALIDATION
   implicit none
   private

   public :: device_memory_free_bytes, device_memory_fits, device_memory_set_cap, device_memory_local_growth

   integer(int64), save :: free_cap = 0_int64  !! the test hook's answer, 0 for none

contains

   subroutine device_memory_set_cap(bytes)
      !! The free device memory the query reports, for tests; 0 or less removes the cap
      integer(int64), intent(in) :: bytes

      free_cap = max(bytes, 0_int64)
   end subroutine device_memory_set_cap

   function device_memory_free_bytes() result(bytes)
      !! The free device memory now; -1 on the host build or when the query gives nothing (no device to fit).
      !! The test hook's cap when one is set.
#ifdef TRC_DEVICE
      use openacc, only: acc_get_property, acc_get_device_num, acc_get_device_type, acc_property_free_memory
#endif
      integer(int64) :: bytes

      bytes = -1
#ifdef TRC_DEVICE
      bytes = int(acc_get_property(acc_get_device_num(acc_get_device_type()), acc_get_device_type(), &
                                   acc_property_free_memory), int64)
      if (bytes <= 0) bytes = -1
#endif
      if (free_cap > 0) bytes = free_cap
   end function device_memory_free_bytes

   function device_memory_local_growth(frame_bytes) result(bytes)
      !! Device memory, not yet taken, that CUDA reserves to launch a kernel whose per-thread stack frame is
      !! `frame_bytes` (read it with `cuobjdump -res-usage` from the linked library): the multiprocessor count
      !! times the resident threads of each, times what the frame exceeds the current stack limit by.
      !! `cudaDeviceGetLimit(cudaLimitStackSize)` reports the grown size after such a launch (measured, NVHPC 26.9,
      !! V100: 1024 before, 53920 after a 53920 B kernel, and the free memory down by 52896 B a thread), and the
      !! free-memory query already counts what is reserved, so only the growth is returned. A host build, or a
      !! failing query, answers as a V100 (80 multiprocessors of 2048 threads) with nothing reserved, so that a
      !! test with `device_memory_set_cap` sees a device's figure.
      integer(int64), intent(in) :: frame_bytes
      integer(int64) :: bytes

#ifdef TRC_DEVICE
      interface
         function cuda_get_device(device) bind(c, name="cudaGetDevice") result(ierr)
            import :: c_int
            implicit none
            integer(c_int), intent(out) :: device
            integer(c_int) :: ierr
         end function cuda_get_device
         function cuda_device_get_attribute(val, attr, device) bind(c, name="cudaDeviceGetAttribute") result(ierr)
            import :: c_int
            implicit none
            integer(c_int), intent(out) :: val
            integer(c_int), value :: attr, device
            integer(c_int) :: ierr
         end function cuda_device_get_attribute
         function cuda_device_get_limit(val, limit) bind(c, name="cudaDeviceGetLimit") result(ierr)
            import :: c_int, c_size_t
            implicit none
            integer(c_size_t), intent(out) :: val
            integer(c_int), value :: limit
            integer(c_int) :: ierr
         end function cuda_device_get_limit
      end interface
      ! cudaDeviceAttr and cudaLimit values of cuda's driver_types.h.
      integer(c_int), parameter :: ATTR_MULTIPROCESSORS = 16, ATTR_THREADS_PER_MULTIPROCESSOR = 39
      integer(c_int), parameter :: LIMIT_STACK_SIZE = 0
      integer(c_int) :: device, nmp, nthread
      integer(c_size_t) :: frame_now
#endif
      integer(int64) :: threads, reserved

      threads = 80_int64*2048_int64
      reserved = 0
#ifdef TRC_DEVICE
      if (cuda_get_device(device) == 0) then
         if (cuda_device_get_attribute(nmp, ATTR_MULTIPROCESSORS, device) == 0) then
            if (cuda_device_get_attribute(nthread, ATTR_THREADS_PER_MULTIPROCESSOR, device) == 0) then
               threads = int(nmp, int64)*int(nthread, int64)
               if (cuda_device_get_limit(frame_now, LIMIT_STACK_SIZE) == 0) reserved = int(frame_now, int64)
            end if
         end if
      end if
#endif
      bytes = max(frame_bytes - reserved, 0_int64)*threads
   end function device_memory_local_growth

   function device_memory_fits(status, what, need, hint) result(fits)
      !! True when `need` bytes fit in the free device memory (always when it is unknown). Otherwise false, and
      !! `status` carries a message that names `what` (start it with the module prefix), the bytes it needs and the
      !! bytes free, in GB, and a remedy: `hint`, one clause, or by default the general one.
      type(error_t), intent(inout) :: status
      character(len=*), intent(in) :: what
      integer(int64), intent(in) :: need
      character(len=*), intent(in), optional :: hint
      logical :: fits

      character(len=*), parameter :: DEFAULT_HINT = "use a smaller basis (or a larger frozen core), or a device with more memory"

      integer(int64) :: free
      character(len=:), allocatable :: remedy

      free = device_memory_free_bytes()
      fits = free < 0 .or. need <= free
      if (fits) return
      remedy = DEFAULT_HINT
      if (present(hint)) remedy = hint
      call status%set(ERROR_VALIDATION, what//" needs "//gb(need)//" GB of device memory and "//gb(free)// &
                      " GB are free: "//remedy)
   end function device_memory_fits

   function gb(bytes) result(text)
      !! `bytes` in GB (1e9), three decimals
      integer(int64), intent(in) :: bytes
      character(len=:), allocatable :: text

      text = to_string_fixed(real(bytes, dp)/1.0e9_dp, 3)
   end function gb

end module trc_device_memory
