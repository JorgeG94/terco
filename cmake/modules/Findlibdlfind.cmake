#
# libdlfind: Kjell Jorner's C-API wrapper around DL-FIND (J. Kaestner, STFC
# Daresbury), fetched with FetchContent the way pic/json-fortran are
# (cmake/modules/sample_utils.cmake's my_fetch_package) -- but NOT through
# that macro directly, and NOT linked into libterco. See TERCO_ENABLE_DLFIND's
# comment in the root CMakeLists.txt for the licensing reason; app/CMakeLists.txt
# is the only place this module is find_package()'d from, and it links the
# result only into terco_cli.
#
# my_fetch_package assumes the fetched project's own CMake target is named
# the same as the package argument ("pic" ships target `pic`, "jsonfortran"
# ships target `jsonfortran`). libdlfind's own CMakeLists.txt does not follow
# that: it names its shared library target `dlfind`, not `libdlfind`
# (repository name and CMake target name differ upstream). Using the macro
# as-is would alias "libdlfind::libdlfind" to a target literally named
# "libdlfind", which does not exist, and the failure would only show up at
# the FINAL LINK of terco_cli, not at configure time. So this Find module
# repeats my_fetch_package's FetchContent_Declare/MakeAvailable shape by
# hand, and aliases to the target libdlfind's CMakeLists actually produces.
#
# Pinned to a COMMIT, not a branch, for the same reason as every other
# fetched dependency here: what terco builds against must not move because
# someone pushed to the fork. 4167998 is upstream's own
# (digital-chemistry-laboratory/libdlfind) HEAD as of 2026-04-06 ("Update
# README.md"); JorgeG94/libdlfind is a fork that exists only so this pin
# cannot move under us. terco carries no patches of its own on top of it
# either.
#
# Offline builds (a compute node with no network): point
# -DFETCHCONTENT_SOURCE_DIR_LIBDLFIND=<a local clone checked out to the
# pinned commit above>.
#
# libdlfind's OWN CMakeLists.txt (not terco's) builds it: a SHARED library
# named `dlfind`, serial (dlf_serial.f90 -- no MPI inside DL-FIND at all),
# linked against BLAS/LAPACK it finds itself, with `-std=legacy` under GNU
# only (DL-FIND's F90 predates gfortran's current strictness; nvfortran
# needs no extra flag). terco adds nothing to that build and declares
# `api_dl_find` itself in app/trc_dlfind_bridge.F90's own bind(C) interface
# block rather than reading libdlfind's .mod files, so the two projects may
# be built by different Fortran compilers.
#
include(FetchContent)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  set(FETCHCONTENT_QUIET FALSE)
endif()

message(STATUS "Retrieving libdlfind from https://github.com/JorgeG94/libdlfind.git "
               "(ref: 4167998d16d8dac4a484ba9305f27d6325a7a28d)")

# EXCLUDE_FROM_ALL: a dependency built for terco_cli to link, not shipped --
# same reasoning as my_fetch_package's own comment (otherwise libdlfind's
# own `install()` rules would install into terco's prefix too).
FetchContent_Declare(
  libdlfind
  GIT_REPOSITORY "https://github.com/JorgeG94/libdlfind.git"
  GIT_TAG "4167998d16d8dac4a484ba9305f27d6325a7a28d"
  EXCLUDE_FROM_ALL)
FetchContent_MakeAvailable(libdlfind)

if(NOT TARGET dlfind)
  message(FATAL_ERROR "libdlfind fetched, but its CMakeLists.txt did not define "
                       "the expected `dlfind` target -- has upstream renamed it?")
endif()

if(NOT TARGET libdlfind::libdlfind)
  add_library(libdlfind::libdlfind INTERFACE IMPORTED)
  target_link_libraries(libdlfind::libdlfind INTERFACE dlfind)
endif()
