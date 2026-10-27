# With TERCO_BUILD_SHARED=OFF, libterco_opt must be static and the terco
# executable must not load a libterco_opt.so. A shared libterco_opt beside
# a static libterco carries its own copy of terco and pic-mpi (and of the
# legacy `use mpi` common blocks), and with Cray MPICH that split
# MPI_IN_PLACE: multi-rank runs hung at the first Fock allreduce. See
# src/opt/CMakeLists.txt and .claude/context/gpu-pitfalls.md.
#
#   cmake -DEXE=<terco> -DOPT=<libterco_opt file> -DOBJDUMP=<objdump> \
#         -P check_static_terco_opt.cmake

if(NOT OPT MATCHES "\\.a$")
  message(FATAL_ERROR "libterco_opt is not a static archive in a TERCO_BUILD_SHARED=OFF build: ${OPT}")
endif()

execute_process(
  COMMAND ${OBJDUMP} -p ${EXE}
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "${OBJDUMP} -p ${EXE} failed (${rc}):\n${err}")
endif()
if(out MATCHES "NEEDED +libterco_opt")
  message(FATAL_ERROR "${EXE} loads libterco_opt.so in a TERCO_BUILD_SHARED=OFF build")
endif()
if(NOT out MATCHES "NEEDED +libdlfind")
  message(FATAL_ERROR "${EXE} does not link libdlfind.so; the check is not looking at the right binary")
endif()
