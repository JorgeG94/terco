#
# toml-f: the TOML parser behind the `terco` CLI's input files (docs/INPUT.md,
# app/trc_run_input.F90). Only terco_cli links it; libterco never does.
#
# Pinned to JorgeG94/toml-f's chore/get_to_work_nvhpc (e0bc8c6), commits on
# upstream f9ff547 that this build needs: -Mbackslash/-Mallocatable=03 under
# NVHPC in toml-f's CMake (upstream does not build through CMake with
# nvfortran), convert_ucs fixed for nvfortran (a four-byte \U escape came back
# padded with blanks), and a cmake_minimum_required policy range (without it
# every configure prints CMake's < 3.10 deprecation warning). Move to an
# upstream tag once those are merged.
#
# Built STATIC and without its tests: toml-f's library type follows
# BUILD_SHARED_LIBS, and an excluded-from-all shared library would not be
# installed next to the terco executable that needs it; its tests follow
# BUILD_TESTING, which would register toml-f's own checks in terco's ctest.
# Both are set only for this fetch and put back afterwards.
#
set(_lib "toml-f")
set(_url "https://github.com/JorgeG94/toml-f.git")
include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")
set(_rev "e0bc8c69985adcae9b61af45b8dd18ceea54ef2c")

set(_terco_saved_shared "${BUILD_SHARED_LIBS}")
set(_terco_saved_testing "${BUILD_TESTING}")
set(BUILD_SHARED_LIBS OFF)
set(BUILD_TESTING OFF)
my_fetch_package("${_lib}" "${_url}" "${_rev}")
set(BUILD_SHARED_LIBS "${_terco_saved_shared}")
set(BUILD_TESTING "${_terco_saved_testing}")

unset(_terco_saved_shared)
unset(_terco_saved_testing)
unset(_lib)
unset(_url)
unset(_rev)
