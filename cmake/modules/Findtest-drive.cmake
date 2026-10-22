#
# test-drive: the unit-test framework new checks are written in (collect a
# suite of named unit tests, `check` with a message, run them and report
# each). Only the test build fetches it; the library never links it.
#
set(_lib "test-drive")
set(_url "https://github.com/fortran-lang/test-drive.git")
include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")
# Three commits past v0.6.1 on main: among them the check for whether the
# compiler has the XDP and QP kinds, without which v0.6.1 does not build
# with Flang. So this pins the exact commit roundabout's CI verified at
# 174/174 on gfortran, ifx, nvfortran and flang. A SHA is a stronger pin than a tag
# anyway (tags can be moved); revisit when test-drive cuts a release newer
# than v0.6.1. (metalquicha pins v0.6.1 -- it does not build with Flang, so
# it never hits this.)
set(_rev "791484a79a24cbc00be6f74e53c6e030cb75efe2")
set(TESTDRIVE_BUILD_TESTING OFF CACHE BOOL "" FORCE)
my_fetch_package("${_lib}" "${_url}" "${_rev}")
unset(_lib)
unset(_url)
unset(_rev)
