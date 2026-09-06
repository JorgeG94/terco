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
if(TERCO_ENABLE_GPU)
  set(TERCO_GPU_OPTS "${TERCO_GPU_ARCH},mem:separate")
  if(TERCO_CUDA_VERSION)
    set(TERCO_GPU_OPTS "${TERCO_GPU_OPTS},cuda${TERCO_CUDA_VERSION}")
  endif()
  set(TERCO_FLAGS "-stdpar=gpu -acc -gpu=${TERCO_GPU_OPTS} -Mpreprocess -traceback")
else()
  #
  # Host build. No `-acc`, so every OpenACC directive is an inert comment, and
  # `do concurrent` is an ordinary loop -- which is exactly why the answers
  # must be unchanged. The atomics are comments too and cannot race, because
  # nothing is running concurrently.
  #
  # Only the preprocessor flag differs by compiler; the sources are otherwise
  # standard Fortran.
  #
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
    # -Mnovect: nvfortran 26.5's vectoriser miscompiles an unrolled loop in
    # the (sd|..) four-centre kernels at -O2 and above ON THE HOST -- NaN in
    # the Fock matrix for any basis with d shells and long contractions
    # (cc-pVDZ water: 169 of them; -O1, -Mnovect or -Mnounroll all clean,
    # gfortran clean at every level). The device build compiles those loops
    # for the GPU and is unaffected, and this host build exists to check
    # numbers, not to be fast.
    set(TERCO_FLAGS "-Mpreprocess -traceback")
    # As a target option, not here: CMake's own Release flags (-fast -O3)
    # come after these and -fast switches the vectoriser back on.
    set(TERCO_HOST_NVHPC_OPTIONS "-Mnovect" PARENT_SCOPE)
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(TERCO_FLAGS "-cpp -ffree-line-length-none -fbacktrace")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel|IntelLLVM")
    set(TERCO_FLAGS "-fpp -traceback")
  elseif(CMAKE_Fortran_COMPILER_ID MATCHES "LLVMFlang|Flang")
    set(TERCO_FLAGS "-cpp")
  else()
    message(
      WARNING
        "Unknown Fortran compiler '${CMAKE_Fortran_COMPILER_ID}'; assuming "
        "-cpp enables the preprocessor. Please add it here.")
    set(TERCO_FLAGS "-cpp")
  endif()
endif()

set(CMAKE_Fortran_FLAGS
    "${CMAKE_Fortran_FLAGS} ${TERCO_FLAGS}"
    PARENT_SCOPE)
