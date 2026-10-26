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
# someone pushed to the fork. c5eb569 is JorgeG94/libdlfind's own
# terco-patches branch, three patches on top of upstream's
# (digital-chemistry-laboratory/libdlfind) 4167998 HEAD (2026-04-06,
# "Update README.md"), that this bridge needs:
#   - 86ca316 passes the `kiter` argument dlf_parallel_opt.f90's three
#     dlf_get_gradient calls were missing (the stochastic/genetic global
#     searches segfaulted before the first energy evaluation without it),
#     and lets DLF_PO_SEED (an integer environment variable) replace the
#     wall-clock seed those searches otherwise reseed from on every call.
#   - 3c27977 allocates glob%xcoords2 and sets tcoords2 before the
#     connect task (1011) stores the P-RFO transition mode into it, which
#     a connect run given no second geometry (every --opt-connect run
#     through this bridge) otherwise leaves unallocated -- the cause of
#     the P-RFO + --opt-connect crash app/main.F90 used to refuse outright.
#   - c5eb569 makes dlf_cartesian_align's eigen-decomposition deterministic
#     for degenerate (linear/planar) alignments. dlf_neb_init Kabsch-aligns
#     every known point against the start geometry by diagonalising a 3x3
#     matrix that is rank-deficient for a linear molecule (HCN/HNC's own
#     NEB endpoints), so the eigensolver was free to return any orthonormal
#     basis of the near-zero-eigenvalue null space -- decided by
#     floating-point noise that measurably shifted under concurrent MKL
#     load. The fix canonicalises the eigenbasis (sign-fixed eigenvectors,
#     and a canonical basis for each cluster of near-equal eigenvalues)
#     right after the diagonalisation, making the alignment a deterministic
#     function of the input geometry; this is what let terco drop the
#     retry-once workaround that used to live in test/check_opt_neb.F90's
#     run_neb.
# GIT_TAG is the commit, not the branch name, so a push to terco-patches
# after this cannot move the build out from under us either; terco-patches
# itself exists only so these three commits have somewhere to live upstream
# of a fork PR. GIT_SHALLOW FALSE below is spelled out explicitly (it is
# already FetchContent_Declare's own default) because it is load-bearing
# here: terco-patches is not libdlfind's default branch, and a shallow
# clone only fetches the tip of whichever ref it resolves, so a full
# clone is what lets the pin above -- reachable from terco-patches, not
# from the default branch -- resolve at all.
#
# Offline builds (a compute node with no network): point
# -DFETCHCONTENT_SOURCE_DIR_LIBDLFIND=<a local clone checked out to the
# pinned commit above>. `git clone` the fork, then `git checkout
# c5eb56901040e831d46d3fc073be19141cc0525f` gets one.
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
               "(ref: c5eb56901040e831d46d3fc073be19141cc0525f, branch terco-patches)")

# EXCLUDE_FROM_ALL: a dependency built for terco_cli to link, not shipped --
# same reasoning as my_fetch_package's own comment (otherwise libdlfind's
# own `install()` rules would install into terco's prefix too).
# GIT_SHALLOW FALSE: see this file's own header for why this is spelled
# out rather than left at its default.
FetchContent_Declare(
  libdlfind
  GIT_REPOSITORY "https://github.com/JorgeG94/libdlfind.git"
  GIT_TAG "c5eb56901040e831d46d3fc073be19141cc0525f"
  GIT_SHALLOW FALSE
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
