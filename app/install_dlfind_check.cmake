# An installed DLFIND=ON terco_cli needs libdlfind.so at its own RPATH
# ($ORIGIN/../lib, app/CMakeLists.txt) -- Findlibdlfind.cmake's
# EXCLUDE_FROM_ALL keeps libdlfind's OWN install() rules from firing, so
# without an explicit install(TARGETS dlfind ...) of our own, the shared
# library is built but never copied into the install prefix, and the
# installed binary fails to even start. This installs to a scratch prefix
# and starts the installed binary for real, rather than only checking a
# file exists:
#
#   cmake -DCMAKE_COMMAND=<cmake> -DBUILD_DIR=<build tree> -DPREFIX=<scratch prefix> \
#         -P install_dlfind_check.cmake

file(REMOVE_RECURSE ${PREFIX})

execute_process(
  COMMAND ${CMAKE_COMMAND} --install ${BUILD_DIR} --prefix ${PREFIX}
  RESULT_VARIABLE rc
  OUTPUT_VARIABLE out
  ERROR_VARIABLE err)
if(NOT rc EQUAL 0)
  message(FATAL_ERROR "cmake --install failed (${rc}):\n${out}\n${err}")
endif()

if(NOT EXISTS ${PREFIX}/bin/terco)
  message(FATAL_ERROR "no ${PREFIX}/bin/terco after install")
endif()

# libdlfind's own CMakeLists.txt builds it as `libdlfind.so.<version>`
# with a `libdlfind.so` symlink (the usual SONAME convention) -- glob
# rather than assume which of the two names app/CMakeLists.txt's
# install(TARGETS dlfind ...) actually copied.
file(GLOB _dlfind_installed ${PREFIX}/lib/libdlfind.so*)
if(NOT _dlfind_installed)
  message(FATAL_ERROR "no libdlfind.so* under ${PREFIX}/lib after install -- "
                       "the dlfind install() rule (app/CMakeLists.txt) did not run")
endif()

# The real test: starting the installed binary away from the build tree,
# with no LD_LIBRARY_PATH help, resolves libdlfind.so through RPATH alone.
execute_process(
  COMMAND ${PREFIX}/bin/terco --version
  RESULT_VARIABLE rc2
  OUTPUT_VARIABLE out2
  ERROR_VARIABLE err2)
if(NOT rc2 EQUAL 0)
  message(FATAL_ERROR "installed terco --version failed (${rc2}), likely could not load "
                       "libdlfind.so through its RPATH:\n${out2}\n${err2}")
endif()

message(STATUS "install check: terco --version ran from ${PREFIX} with libdlfind resolved "
               "through $ORIGIN/../lib")
