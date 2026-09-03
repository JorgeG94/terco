cmake_minimum_required(VERSION 3.28)

project(
  "TERCO"
  LANGUAGES Fortran
  VERSION 0.1.0
  DESCRIPTION "GPU Gaussian integrals in standard Fortran")
set(project_name terco)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
include(GNUInstallDirs)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/modules")
# Named for the project, not `print_status`: PIC ships a module by that name
# and does `include(print_status)` itself. With terco's module path searched
# first, a shared name means PIC includes OURS and then cannot find its own
# function.
include(terco_print_status)

# ---------------------------------------------------------------------------
# Options
# ---------------------------------------------------------------------------

#
# The four-centre angular momentum ceiling.
#
# Not a preference. The four-centre kernels are generated one per (la,lb,lc,ld)
# class, so the count goes as L^4: LMAX=1 is 16 classes, LMAX=2 is 81 and
# (dd|dd) alone is twenty thousand lines. f would need both a larger Boys table
# -- BOYS_MMAX caps la+lb+lc+ld at 8 and (ff|ff) needs 12 -- and an answer to
# compile time.
#
# 2 is the shipped configuration and the one whose generated sources are in the
# repository. Anything else is regenerated into the build tree, which needs
# Python 3.
#
set(TERCO_LMAX
    "2"
    CACHE STRING "Four-centre angular momentum ceiling (0=s, 1=p, 2=d)")
set_property(CACHE TERCO_LMAX PROPERTY STRINGS 0 1 2)

set(TERCO_GPU_ARCH
    "cc70"
    CACHE STRING "Target GPU architecture, as nvfortran spells it (cc70, cc80, ...)")

option(TERCO_ENABLE_TESTING "Build the correctness checks" ON)
option(TERCO_REGENERATE
       "Regenerate the committed kernels in-tree rather than using them" OFF)
#
# Shared by default. The one consumer that matters, metalquicha, is built
# with gfortran and reaches terco through the C ABI, and a static libterco
# from nvfortran cannot be linked by another compiler: its device code, the
# OpenACC runtime and now cuBLAS and cuSOLVER all expect nvfortran's own
# link step. A shared library takes that step once, here, with the rpath
# baked in, and hands the consumer an ordinary .so with a C symbol table.
# OFF is for terco's own tests and nvfortran consumers, where either works.
#
option(TERCO_BUILD_SHARED
       "Build a shared library, for host codes built with another compiler" ON)

#
# The GPU build needs NVHPC: `-stdpar=gpu` is what turns `do concurrent` into
# device kernels, and no other compiler does that today.
#
# A CPU build is still worth having, and not as a fallback anyone should run
# production on. It is the only configuration a CI runner without a GPU can
# compile and VALIDATE, and validating on the host is what makes a green CI
# mean something: `do concurrent` is serial there, the OpenACC directives are
# inert comments, and the integrals must come out the same. A number that
# depends on the offload model is a bug in the offload model.
#
option(TERCO_ENABLE_GPU "Offload with -stdpar=gpu (needs NVHPC)" ON)

if(TERCO_ENABLE_GPU AND NOT CMAKE_Fortran_COMPILER_ID STREQUAL "NVHPC")
  message(
    FATAL_ERROR
      "TERCO_ENABLE_GPU=ON requires the NVHPC Fortran compiler (found "
      "'${CMAKE_Fortran_COMPILER_ID}'). Configure with "
      "-DCMAKE_Fortran_COMPILER=nvfortran, or -DTERCO_ENABLE_GPU=OFF for a "
      "host build.")
endif()

add_subdirectory(cmake)

#
# BEFORE anything creates a target, PIC included.
#
# CMake reads this when a target is CREATED, to initialise that target's
# POSITION_INDEPENDENT_CODE property -- so setting it after `find_package(pic)`
# leaves PIC's static library built without -fPIC, and linking it into a shared
# libterco.so then fails:
#
#     relocation R_X86_64_32 against `.rodata' can not be used when making a
#     shared object; recompile with -fPIC
#
# Which only shows up with TERCO_BUILD_SHARED=ON, i.e. exactly the
# configuration a host code built with another compiler needs.
#
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#
# PIC, for one thing: `sort_index` on 64-bit Morton codes, to order the shell
# pairs inside each bin spatially. See `sort_bins_spatially` in
# `src/trc_bins.F90` for why that ordering is worth a sort.
#
# PIC carries these sorts precisely because the Fortran stdlib does not build
# with the NVIDIA compilers, which is the situation terco is permanently in.
#
#
# PIC's own tests are not terco's business. Fetched with them on, `ctest` here
# runs them too -- and three of them fail under `-fopenmp` on the host for
# reasons that have nothing to do with integrals, which would leave this
# project's CI red over a dependency's logger.
#
set(PIC_ENABLE_TESTING
    OFF
    CACHE BOOL "" FORCE)
find_package("pic" REQUIRED)

set(main_lib terco)
if(TERCO_BUILD_SHARED)
  add_library(${main_lib} SHARED)
else()
  add_library(${main_lib} STATIC)
endif()
add_library(terco::${project_name} ALIAS ${main_lib})

#
# TRC_FOCK6 folds the eight-fold permutational symmetry into six atomic
# updates per quartet; TRC_PERCLASS launches one kernel per angular-momentum
# class so ptxas budgets registers for that class rather than the worst one.
# Both are measured wins and neither is optional in practice, so they are
# definitions rather than options -- an option nobody should turn off is a
# configuration that nobody tests.
#
target_compile_definitions(${main_lib} PRIVATE TRC_LMAX=${TERCO_LMAX}
                                               TRC_FOCK6 TRC_PERCLASS)

#
# The unrolled VRR/HRR go into the SHARED kernel only, and only at LMAX<=1.
#
# That kernel holds two copies of the recurrence, so unrolling it at LMAX=2
# puts twice the entire per-class kernel set into one translation unit and
# nvfortran does not finish. Under TRC_PERCLASS the shared kernel carries only
# the non-symmetric Fock build, while every SCF quartet goes through a
# per-class kernel that is unaffected -- so the 2.4x the unroll is worth there
# buys nothing an SCF can see.
#
if(TERCO_LMAX LESS_EQUAL 1)
  target_compile_definitions(${main_lib} PRIVATE TRC_UNROLL_VRR TRC_UNROLL_HRR)
endif()

add_subdirectory(src)

#
# cuSOLVER and cuBLAS, for the SCF's eigenproblem and products
# (src/trc_linalg.F90). A LINK option on terco's own target, PUBLIC so that
# anything linking libterco with nvfortran gets it too -- and NOT in
# CMAKE_Fortran_FLAGS, where it reached every target in the tree: pic's own
# `app`, fetched as a dependency, was then asked to link CUDA math libraries
# it has no use for, and failed wherever the NVHPC install lacks them for
# the driver's CUDA version. The NVHPC modules the source `use`s come from
# the compiler's include path and need nothing at compile time.
#
if(TERCO_ENABLE_GPU)
  target_link_options(${main_lib} PUBLIC -cudalib=cublas,cusolver)
endif()

#
# pic-blas, for the SCF driver (src/trc_scf.F90): the orthogonaliser and the
# core guess go through its explicit-interface wrappers in every build, and
# under a host build so do the eigenproblem, the products and the dot
# products of the iteration. The GPU build does those on the device;
# -cudalib is a link option on the target, above. The DIIS system is solved in terco itself, the way
# metalquicha does it, and needs neither.
#
# pic-blas's switch is the unprefixed ENABLE_TESTING; terco's own is
# TERCO_ENABLE_TESTING, so forcing this off costs nothing here.
set(ENABLE_TESTING
    OFF
    CACHE BOOL "" FORCE)
find_package("pic-blas" REQUIRED)
# Build-interface only. pic and pic-blas are static and end up inside
# libterco.so; an installed consumer links nothing else and the export set
# must not name targets that are not installed. terco's own tests, which
# `use pic_types`, get them through this build interface.
target_link_libraries(${main_lib} PUBLIC $<BUILD_INTERFACE:pic::pic>
                                         $<BUILD_INTERFACE:pic-blas::pic-blas>)

set_target_properties(${main_lib} PROPERTIES Fortran_MODULE_DIRECTORY
                                             ${CMAKE_BINARY_DIR}/modules)
target_include_directories(
  ${main_lib} PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
                     $<INSTALL_INTERFACE:include/terco>)

if(TERCO_ENABLE_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

# ---------------------------------------------------------------------------
# Install and export
# ---------------------------------------------------------------------------

install(
  TARGETS ${main_lib}
  EXPORT tercoTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES
  DESTINATION include/terco)

# terco's own modules only. The modules/ directory is shared with the
# fetched dependencies, whose .mod files are theirs to ship; and the test
# helper's module is not part of the library.
install(
  DIRECTORY ${CMAKE_BINARY_DIR}/modules/
  DESTINATION include/terco
  FILES_MATCHING
  PATTERN "trc_test_basis.mod" EXCLUDE
  PATTERN "trc_*.mod"
  PATTERN "bit_repro*.mod")

#
# The C ABI's Fortran declaration ships as SOURCE, not as a module.
#
# A consumer built with another compiler cannot read terco's `.mod` files, so
# it compiles this one file itself and gets the interfaces from the same
# declaration terco's own tests use. Shipping a `.mod` instead would work for
# NVHPC callers only, which are the ones that do not need it.
#
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/include/trc_c_interfaces.f90
        DESTINATION include/terco)

install(
  EXPORT tercoTargets
  FILE tercoTargets.cmake
  NAMESPACE terco::
  DESTINATION lib/cmake/terco)

include(CMakePackageConfigHelpers)

configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/tercoConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/tercoConfig.cmake"
  INSTALL_DESTINATION lib/cmake/terco)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/tercoConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/tercoConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/tercoConfigVersion.cmake"
        DESTINATION lib/cmake/terco)

print_terco_summary()
