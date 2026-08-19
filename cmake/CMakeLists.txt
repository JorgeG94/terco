if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Build type to be used." FORCE)
  message(
    STATUS "Setting build type to '${CMAKE_BUILD_TYPE}' as none was specified.")
endif()

#
# The GPU flags, which are the build.
#
#   -stdpar=gpu   `do concurrent` becomes device kernels. This is the whole
#                 mechanism; without it the same source is a slow CPU library.
#   -acc          the OpenACC data directives and atomics that `do concurrent`
#                 cannot express -- residency, and the scatter into six
#                 overlapping blocks per quartet.
#   mem:separate  NOT managed memory. Host and device hold distinct copies and
#                 every array a kernel touches is mapped explicitly, so a
#                 missing map is a link-time or run-time failure rather than a
#                 silent page fault that costs an order of magnitude.
#   -Mpreprocess  the sources are .F90 and select on TRC_LMAX.
#
set(TERCO_GPU_FLAGS
    "-stdpar=gpu -acc -gpu=${TERCO_GPU_ARCH},mem:separate -Mpreprocess -traceback")

set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} ${TERCO_GPU_FLAGS}"
    PARENT_SCOPE)
