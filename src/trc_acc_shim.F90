!
! A no-op `acc_register_library`, so terco can be loaded beside GNU libgomp.
!
! THE BUG THIS FIXES
! ------------------
! NVHPC's OpenACC runtime, `libacchost.so`, carries a **weak undefined**
! reference to `acc_register_library` -- the OpenACC profiling interface's
! registration hook. Weak and undefined means: if some profiling tool defines
! it, call it; if nothing does, it resolves to null and the runtime skips it.
! In an all-NVHPC process nothing defines it, and that is the normal path.
!
! GNU `libgomp` DEFINES it, as an unimplemented stub that prints
!
!     libgomp: TODO
!
! and kills the process.
!
! So the moment terco is linked into a program that also uses GNU OpenMP --
! which is any gfortran host code built with `-fopenmp`, mqc among them --
! libacchost's weak reference binds to libgomp's stub and the program dies at
! the first OpenACC call, before any integral is computed. Nothing in the
! message names OpenACC, terco, or NVHPC.
!
! WHY THE FIX IS HERE AND WHY IT IS A NO-OP
! ------------------------------------------
! Defining the symbol in `libterco.so` wins, because the host links terco
! directly and the dynamic linker takes the first definition in DT_NEEDED
! order -- terco is a direct dependency, libgomp's definition is found later.
!
! A no-op is not a workaround that papers over something: it is EXACTLY the
! behaviour of the weak symbol resolving to null, which is what happens
! whenever libgomp is not in the process. The only thing given up is the
! ability of a third-party OpenACC profiling tool to register itself through
! libgomp, which was never going to work anyway -- libgomp's version of it
! aborts.
!
! This belongs to terco rather than to any host code because terco is what
! brings the NVHPC OpenACC runtime into a process it did not previously
! inhabit. A host code cannot reasonably be expected to know that linking a
! GPU library changes how its OpenMP runtime is consulted.
!
module trc_acc_shim
   use, intrinsic :: iso_c_binding
   implicit none
   private
   public :: acc_register_library

contains

   !> The OpenACC profiling registration hook, deliberately doing nothing.
   !>
   !> The three arguments are the register, unregister and lookup callbacks
   !> the OpenACC profiling interface passes; they are taken by value as
   !> opaque pointers and ignored.
   subroutine acc_register_library(reg, unreg, lookup) &
      bind(c, name="acc_register_library")
      type(c_ptr), value :: reg, unreg, lookup
      ! Reference them so no compiler warns about unused dummies while
      ! keeping the body genuinely empty.
      if (c_associated(reg) .and. c_associated(unreg) &
          .and. c_associated(lookup)) return
   end subroutine acc_register_library

end module trc_acc_shim
