#
# libterco_opt: terco's geometry optimiser (DL-FIND, via libdlfind), as its
# own shared library -- see the root CMakeLists.txt's TERCO_ENABLE_DLFIND
# option for the licensing reason it is not part of ${main_lib}.
#
# This directory is only added (root CMakeLists.txt) when
# TERCO_ENABLE_DLFIND is ON. With it OFF, none of this runs: no libdlfind
# fetch, no terco_opt target, nothing installed -- app/CMakeLists.txt's own
# trc_dlfind_stub.F90 is the whole --opt story in that build.
#
# WHY A SEPARATE LIBRARY AND NOT PART OF ${main_lib}. Two reasons, same as
# before this was split out of app/trc_dlfind_bridge.F90:
#   - libdlfind is LGPL-3.0-or-later; libterco is not. Keeping DL-FIND
#     out of libterco.so is what keeps libterco's own licence clean and
#     the LGPL relink obligation scoped to this one library.
#   - DL-FIND's own `error stop` (api.f90's dlf_error) can end the WHOLE
#     PROCESS, never just report a status code back -- a library that
#     might do that must be something a caller opts into explicitly
#     (TERCO_ENABLE_DLFIND=ON, and then linking terco::terco_opt), never
#     something pulled in by linking libterco for an energy or a gradient.
#
# WHO LINKS THIS. terco's own CLI (app/CMakeLists.txt, below this
# directory in the build) and cyclone or any other out-of-tree C/C++/
# Fortran caller, through find_package(terco) and terco::terco_opt
# (the install/export block at the bottom) plus include/terco_opt.h.
# Both reach the exact same Fortran module `trc_dlfind_bridge` --
# terco's own CLI calls trc_dlfind_optimize directly (app/main.F90, a
# Fortran caller built by the same compiler); a C caller goes through
# trc_opt_capi.F90's bind(c) skin instead. One driver, not a copy per
# caller.
#
find_package(libdlfind REQUIRED)

add_library(terco_opt SHARED trc_opt_io.F90 trc_dlfind_bridge.F90 trc_opt_capi.F90)
add_library(terco::terco_opt ALIAS terco_opt)

# libterco_opt.so, not terco_cli, now carries the DT_NEEDED on
# libdlfind.so, and an executable's RUNPATH does not cover its libraries'
# own dependencies, so this library needs its own: libdlfind.so is
# installed next to it in lib/.
if(APPLE)
  set_target_properties(terco_opt PROPERTIES INSTALL_RPATH "@loader_path")
else()
  set_target_properties(terco_opt PROPERTIES INSTALL_RPATH "$ORIGIN")
endif()

set_target_properties(terco_opt PROPERTIES Fortran_MODULE_DIRECTORY
                                           ${CMAKE_BINARY_DIR}/modules)
target_include_directories(
  terco_opt PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules>
                   $<INSTALL_INTERFACE:include/terco>)

# PUBLIC main_lib: trc_dlfind_bridge.F90/trc_opt_capi.F90 `use` trc_capi
# (part of ${main_lib}), so a consumer linking terco::terco_opt needs
# terco's own .mod search path and symbols too, not just this library's --
# PUBLIC on this link line is what makes `target_link_libraries(x PRIVATE
# terco::terco_opt)` pull terco::terco in for `x` automatically.
#
# PRIVATE libdlfind::libdlfind: no libdlfind type or symbol crosses this
# library's own public API (trc_opt_capi.F90's handles are opaque, and
# api_dl_find's interface is declared by hand in trc_dlfind_bridge.F90,
# never through libdlfind's own .mod -- see that file's own header). A
# shared library's PRIVATE link requirement still becomes a DT_NEEDED
# entry on libterco_opt.so (the symbols ARE used, just not re-exposed), so
# a consumer still needs libdlfind.so found at RUNTIME -- never at
# configure/link time, which is the whole point of not exposing it.
target_link_libraries(terco_opt PUBLIC ${main_lib} PRIVATE libdlfind::libdlfind pic::pic)

install(
  TARGETS terco_opt
  EXPORT tercoTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES
  DESTINATION include/terco)

# libdlfind.so itself, alongside libterco_opt.so and libterco.so -- an
# installed TERCO_ENABLE_DLFIND=ON terco could not start at all without
# this (Findlibdlfind.cmake's own EXCLUDE_FROM_ALL keeps libdlfind's OWN
# install() rules, upstream, from firing on its own). terco_cli's RPATH
# ($ORIGIN/../lib, app/CMakeLists.txt) already covers this lib/, same as
# libterco_opt.so's own.
install(TARGETS dlfind LIBRARY DESTINATION lib RUNTIME DESTINATION lib)

# terco_opt's own modules, same convention as the root CMakeLists.txt's
# install of terco's.
install(
  DIRECTORY ${CMAKE_BINARY_DIR}/modules/
  DESTINATION include/terco
  FILES_MATCHING
  PATTERN "trc_dlfind_bridge.mod"
  PATTERN "trc_opt_io.mod"
  PATTERN "trc_opt_capi.mod")

# The C ABI's Fortran declaration ships as SOURCE, not as a module, and its
# C header ships alongside it -- same reasoning as the root CMakeLists.txt's
# install of trc_c_interfaces.f90/terco.h.
install(FILES ${PROJECT_SOURCE_DIR}/include/trc_opt_c_interfaces.f90
              ${PROJECT_SOURCE_DIR}/include/terco_opt.h
        DESTINATION include/terco)
