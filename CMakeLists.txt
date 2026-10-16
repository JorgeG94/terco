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
# (dd|dd) alone is twenty thousand lines. BOYS_MMAX caps la+lb+lc+ld at 16,
# which covers g (gg|gg needs 16), so f and g no longer need a larger Boys
# table -- the remaining obstacle past 2 is compile time.
#
# HGP stops at d. LMAX=3 keeps HGP's d set exactly as at LMAX=2 and builds
# every class with an f shell as a Rys kernel only (trc_rys_select.F90's
# HGP_LMAX).
#
# 2 is the shipped configuration and the one whose generated sources are in the
# repository. Anything else is regenerated into the build tree, which needs
# Python 3.
#
set(TERCO_LMAX
    "2"
    CACHE STRING "Four-centre angular momentum ceiling (0=s, 1=p, 2=d, 3=f)")
set_property(CACHE TERCO_LMAX PROPERTY STRINGS 0 1 2 3)

#
# The density-fitting kernels' orbital ceiling (gen_df.py --lmax; the
# auxiliary ceiling is one higher). Empty follows TERCO_LMAX. Separate
# because DF is three-centre with an auxiliary shell one above the orbital
# one: gen_df.py now emits one module per class (src/CMakeLists.txt's
# --split, same reason gen_perclass.py's does), but before that split, at
# TERCO_LMAX=3 (aux g) nvfortran spent over two hours on the single
# monolithic file and had not finished, while nothing in the four-centre f
# path uses it -- raising this alone (src/CMakeLists.txt sends it to the
# build tree, never src/generated) is how the per-class split is exercised
# without pulling in TERCO_LMAX's unrelated Rys four-centre sweep. 2 keeps
# the shipped DF set; RI-MP2 and the trc_df_* entry points then stop on an
# orbital basis above d.
#
set(TERCO_DF_LMAX
    ""
    CACHE STRING "Density-fitting orbital angular momentum ceiling (empty: TERCO_LMAX)")

#
# TESTING ONLY -- never set by a shipped build. gen_df.py routes a class
# above gen_df.py's SHIPPED_L/SHIPPED_LA (la > 2, lb > 2 or lc > 3 for a
# three-centre class; la > 3 or lc > 3 for a two-centre one) to a Rys
# quadrature body instead of the written-out VRR+HRR -- see that script's
# module docstring. This forces the Rys body onto classes that would
# otherwise keep the VRR+HRR, so the two can be compared element-wise on
# exactly the same (shell pair, auxiliary shell): 'all' forces every class,
# or a comma-separated list of gen_df.py's own class-digit tags (--force-rys
# --help). Threaded straight through to gen_df.py --force-rys; empty means
# "don't pass the flag at all", which is what every normal configuration
# does.
#
set(TERCO_DF_FORCE_RYS
    ""
    CACHE STRING "TESTING ONLY: force gen_df.py's Rys body onto these classes ('' / 'all' / tag list)")
mark_as_advanced(TERCO_DF_FORCE_RYS)

#
# The within-shell Cartesian component order every AO-indexed quantity terco
# produces uses: Fock/J/K, the one-electron matrices, multipoles, DF tensors,
# AO values on the grid, and the density blocks read back in. ORDER ONLY --
# this does not apply GAMESS's per-component sqrt3/sqrt5 normalisation, which
# is a separate decision at the API boundary.
#
# lexicographic (the default) is libcint's own order -- the exponent triples
# in descending lexicographic order, x < y < z alphabetically on the labels
# -- and it is exactly today's behaviour: a lexicographic build needs no new
# -D flag anywhere (src/CMakeLists.txt), so the shipped src/generated and an
# `fpm build`, which passes none, are unaffected. gamess switches d and f
# (s and p are unchanged, being identical in both orders) to GAMESS's own
# order; see scripts/cartorder.py for both tables and their sourcing.
#
# "libcint" is accepted as an alias, same as cartorder.py's own ORDERS dict
# (src/CMakeLists.txt normalises it), so a command line written against the
# old name still works.
#
set(TERCO_CART_ORDER
    "lexicographic"
    CACHE STRING "Within-shell Cartesian component order (lexicographic, gamess)")
set_property(CACHE TERCO_CART_ORDER PROPERTY STRINGS lexicographic gamess libcint)

set(TERCO_GPU_ARCH
    "cc70"
    CACHE STRING "Target GPU architecture, as nvfortran spells it (cc70, cc80, ...)")

#
# The CUDA version the device code is built for, and with it the directory
# `-cudalib` takes cuBLAS and cuSOLVER from: nvfortran resolves
# `-cudalib=cublas,cusolver` to `$NVHPC_ROOT/math_libs/<version>/lib64`, and
# when its default version -- the driver's, or NVHPC_CUDA_HOME's -- has no
# such directory, the link fails with "cannot find -lcusolver" and no path
# named. Perlmutter's NVHPC 26.5 defaults to CUDA 13.2 and ships the math
# libraries for other versions. Empty leaves the compiler's default; a
# version, spelled as `ls $NVHPC_ROOT/math_libs` shows it (12.4, 12.9),
# becomes `-gpu=...,cudaX.Y`, which selects code generation and the math
# libraries together, so the two cannot disagree.
#
set(TERCO_CUDA_VERSION
    ""
    CACHE STRING "CUDA version to build against (empty: nvfortran's default), e.g. 12.9")

option(TERCO_ENABLE_TESTING "Build the correctness checks" ON)
option(TERCO_BUILD_APP "Build and install the terco command line program" ON)
#
# Ranks. The Fock build splits the sorted quartet list across ranks -- every
# nranks-th item from the rank -- and sums the result; one device per rank.
# OFF builds pic-mpi's single-rank backend behind the same API, so the code
# is identical either way and a build without MPI is rank 0 of 1.
#
option(TERCO_ENABLE_MPI "Split the Fock build over MPI ranks, one device each" OFF)
#
# Two different "regenerate" switches:
#
# TERCO_REGENERATE rewrites the COMMITTED kernels in src/generated in place,
# from the generators, so that a diff shows whether the committed set still
# reproduces. It always runs its generators.
#
# TERCO_FORCE_REGENERATE is about the BUILD tree. Every configure
# regenerates what it needs there (the Rys kernels; the whole kernel set
# when TERCO_LMAX, TERCO_DF_LMAX, TERCO_CART_ORDER or TERCO_DF_FORCE_RYS
# leave the shipped configuration) only when a generation stamp says the
# generators, their inputs or their options changed (src/CMakeLists.txt,
# GENERATION STAMPS). ON ignores the stamps for one configure and then
# turns itself back OFF. Editing a generator needs neither: the next build
# reconfigures and regenerates on its own.
#
option(TERCO_REGENERATE
       "Regenerate the committed kernels in-tree rather than using them" OFF)
option(TERCO_FORCE_REGENERATE
       "Regenerate the build-tree kernels once, ignoring the generation stamps" OFF)
#
# Rys quadrature is always built: it is the only source of the analytic
# two-electron gradient (eri%gradient) and every class above d. HGP stays
# the energy path for every class it has unless TRC_RYS_CLASSES routes one
# to Rys at runtime (trc_rys_select.F90). Its cost is compile time --
# roughly 25 CPU-minutes of extra nvfortran GPU compilation, measured per
# object (the chunked heavy gradient classes run to 1-3 min each) -- and
# Python 3 at configure time, which generates the per-class kernels into
# the build tree.
#
# Unroll the Rys kernels' root and component loops (scripts/rys_inline.py's
# UNROLLED MODE) for every class whose Cartesian component count na*nb*nc*nd
# is at most this -- passed straight through to gen_rys.py/gen_rys_grad.py's
# --unroll-max-comp. A class above the threshold keeps today's array/loop
# form. "0" (not unrolling anything) is only useful as a measurement
# baseline; "inf"/"all" unrolls every class, including (dd|dd), which is a
# compile-time measurement point, not a build anyone should ship.
#
# 81 is MEASURED, not a guess (w172, 6-31G*, one GPU): the gradient drops
# from 14.74 s to 10.76 s at 1e-10 (-27%) and 30.59 s to 22.67 s at 1e-13
# (-26%) going from 0 to 81, every class below the threshold matching its
# own e31f06a timing exactly where it is NOT yet unrolled (a built-in
# regression check). Raising the ceiling to 243 unrolls three more classes
# ((dp|pp), (dp|ds), (ds|dd)) but changes the measured total NOT AT ALL
# (10.91 s / 22.89 s, within noise of 81's number) while costing 330-530 s
# of EXTRA serial compile time PER CLASS for those three alone -- compile
# cost that buys nothing. 81 is the point past which unrolling stops paying
# for itself.
#
set(TERCO_RYS_UNROLL_MAX
    "81"
    CACHE STRING "Unroll Rys classes with na*nb*nc*nd <= this ('inf'/'all': every class)")
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
#
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
# `src/eri/trc_bins.F90` for why that ordering is worth a sort.
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
# TRC_CART_ORDER: 0 = lexicographic, 1 = gamess. Defined ONLY for gamess, left
# UNDEFINED for the default: every guard this feature adds reads as
# `#if TRC_CART_ORDER != 0`, and the C preprocessor treats an undefined
# identifier in an `#if` as 0, so a lexicographic build -- fpm's included,
# which passes no -D flags at all -- needs no new define to pass every guard.
# The "libcint" alias is normalised here too (src/CMakeLists.txt does its own
# copy of this same normalisation, for the generator calls it makes below
# `add_subdirectory(src)`; both read the same TERCO_CART_ORDER cache entry).
#
if(TERCO_CART_ORDER STREQUAL "libcint")
  set(TERCO_CART_ORDER "lexicographic")
endif()
if(TERCO_CART_ORDER STREQUAL "gamess")
  target_compile_definitions(${main_lib} PRIVATE TRC_CART_ORDER=1)
elseif(NOT TERCO_CART_ORDER STREQUAL "lexicographic")
  message(FATAL_ERROR "TERCO_CART_ORDER must be lexicographic or gamess, "
                      "got '${TERCO_CART_ORDER}'")
endif()
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
# (src/core/trc_linalg.F90). A LINK option on terco's own target, PUBLIC so that
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
# pic-blas, for the SCF driver (src/scf/trc_scf_driver.F90): the orthogonaliser and the
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
set(PIC_ENABLE_MPI
    ${TERCO_ENABLE_MPI}
    CACHE BOOL "" FORCE)
# Not forced: mpi_f08 is the default, but some MPIs need the legacy `mpi`
# module (Cray MPICH on Perlmutter: -DPIC_USE_LEGACY_MPI=ON
# -DPIC_MPICH_PERLMUTTER=ON), and the -D on the command line has to win.
set(PIC_USE_LEGACY_MPI
    OFF
    CACHE BOOL "Use pic-mpi's legacy mpi-module backend instead of mpi_f08")
find_package("pic-mpi" REQUIRED)
find_package("jsonfortran" REQUIRED)
if(DEFINED TERCO_HOST_NVHPC_OPTIONS)
  target_compile_options(${main_lib} PRIVATE ${TERCO_HOST_NVHPC_OPTIONS})
endif()
#
# The library itself no longer branches on whether MPI is present: every
# pic-mpi backend (mpi_f08, legacy MPI, and the single-rank serial stand-in
# PIC_ENABLE_MPI=OFF selects above) exports the same comm_from_handle, so
# trc_capi's communicator handling compiles and runs unchanged either way.
# The tests do not branch either: they get the Fortran handle they hand the
# C API from pic-mpi's comm_t%fortran_handle(), never from `use mpi_f08`,
# which some MPIs (Cray MPICH under nvfortran) do not provide.
#
target_link_libraries(${main_lib} PUBLIC $<BUILD_INTERFACE:pic::pic>
                                         $<BUILD_INTERFACE:pic-blas::pic-blas>
                                         $<BUILD_INTERFACE:pic-mpi::pic-mpi>
                                         $<BUILD_INTERFACE:jsonfortran::jsonfortran>)

set_target_properties(${main_lib} PROPERTIES Fortran_MODULE_DIRECTORY
                                             ${CMAKE_BINARY_DIR}/modules)
target_include_directories(
  ${main_lib} PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
                     $<INSTALL_INTERFACE:include/terco>)

if(TERCO_ENABLE_TESTING)
  enable_testing()
  add_subdirectory(test)
endif()

# The `terco` command line program (app/), installed alongside the library.
if(TERCO_BUILD_APP)
  add_subdirectory(app)
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
