macro("my_fetch_package" package url)
  # Optional third argument for git reference (branch, tag, or commit)
  set(_git_ref "${ARGV2}")
  if(NOT _git_ref)
    set(_git_ref "HEAD")
  endif()
  # Optional fourth and fifth arguments: a patch file (`patch -p1`) applied
  # to the fetched source right after it is cloned -- for a fix that has
  # not landed upstream yet -- and a symbol it adds, to tell an
  # already-patched checkout from one that isn't on a reconfigure. See
  # Findpic-mpi.cmake's use of this for why one exists.
  set(_patch_file "${ARGV3}")
  set(_patch_marker "${ARGV4}")

  string(TOLOWER "${package}" _pkg_lc)
  string(TOUPPER "${package}" _pkg_uc)

  message(STATUS "Retrieving ${package} from ${url} (ref: ${_git_ref})")
  include(FetchContent)

  if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(FETCHCONTENT_QUIET FALSE)
  endif()

  # EXCLUDE_FROM_ALL: a dependency is built for terco to link, not shipped.
  # Without it every fetched project installs into terco's prefix -- its
  # library, its demo executable, its CMake helpers, and the whole shared
  # modules/ directory once more under its own include/<name>.
  if(_patch_file)
    if(NOT _patch_marker)
      message(FATAL_ERROR "my_fetch_package(${package} ... ${_patch_file}): "
                          "needs a 5th argument naming a symbol the patch adds, "
                          "to tell an already-patched checkout from one that isn't")
    endif()
    find_program(_patch_exe patch REQUIRED)
    find_program(_sh_exe sh REQUIRED)
    # Idempotent: a reconfigure re-runs FetchContent's update step even when
    # it changes nothing (GIT_TAG pinned to a commit), which re-triggers
    # PATCH_COMMAND on a source tree the FIRST configure already patched --
    # and `patch` refuses a hunk it finds already applied. Grepping for a
    # symbol the patch adds is the check: present, nothing to do; absent,
    # apply for real and let a genuine failure stop the build.
    FetchContent_Declare(
      "${_pkg_lc}"
      GIT_REPOSITORY "${url}"
      GIT_TAG "${_git_ref}"
      PATCH_COMMAND "${_sh_exe}" -c
                    "grep -qrF '${_patch_marker}' src/ || ${_patch_exe} -p1 -i ${_patch_file}"
      EXCLUDE_FROM_ALL)
    unset(_patch_exe)
    unset(_sh_exe)
  else()
    FetchContent_Declare(
      "${_pkg_lc}"
      GIT_REPOSITORY "${url}"
      GIT_TAG "${_git_ref}"
      EXCLUDE_FROM_ALL)
  endif()
  FetchContent_MakeAvailable("${_pkg_lc}")

  # Only create the imported interface if the target doesn't already exist This
  # allows the fetched package to create its own namespaced target
  if(NOT TARGET "${package}::${package}")
    add_library("${package}::${package}" INTERFACE IMPORTED)
    target_link_libraries("${package}::${package}" INTERFACE "${package}")
  endif()

  if(NOT EXISTS "${${_pkg_lc}_BINARY_DIR}/include")
    file(MAKE_DIRECTORY "${${_pkg_lc}_BINARY_DIR}/include")
  endif()

  unset(_pkg_lc)
  unset(_pkg_uc)
  unset(_git_ref)
  unset(_patch_file)
  unset(_patch_marker)

  # sanity check
  if(NOT TARGET "${package}::${package}")
    message(FATAL_ERROR "Could not find dependency ${package}")
  endif()
endmacro()
