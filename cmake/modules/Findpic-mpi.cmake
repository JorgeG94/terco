#
# pic-mpi: one API over mpi_f08, the legacy MPI module and a single-rank
# serial backend, so terco depends on it unconditionally and TERCO_ENABLE_MPI
# only chooses the backend (through PIC_ENABLE_MPI). No #ifdef in terco.
#
# PIC_MPI_GIT_TAG is a cache variable rather than a literal below so the
# pinned revision can be switched -- to pic-mpi's comm_from_handle branch
# while it is under review, back to a release tag once it lands on main --
# from the CMake command line (-DPIC_MPI_GIT_TAG=...) without editing this
# file. It currently points at the commit that gives terco
# `comm_from_handle`, needed to accept a caller's own communicator instead
# of only MPI_COMM_WORLD; move it to a tag on main once that branch merges.
#
# Pinned to a COMMIT, not the branch name, so the revision cannot move under
# a build. 3374c5b is pic-mpi's main as of 2b7c247 (0.6.3, which pins pic
# v0.9.2) plus comm_from_handle -- so this is also a step up from the v0.6.0
# terco pinned before, not just one added function. The module-directory
# fix below is the one consequence of that step terco had to absorb.
#
set(PIC_MPI_GIT_TAG
    "3374c5bf58c21c1fd744bf97a570b52882963687"
    CACHE STRING "pic-mpi git ref (branch, tag or commit) to fetch")
set(_lib "pic-mpi")
set(_url "https://github.com/JorgeG94/pic-mpi/")
include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")
my_fetch_package("${_lib}" "${_url}" "${PIC_MPI_GIT_TAG}")
unset(_lib)
unset(_url)

#
# pic-mpi deliberately writes its OWN .mod files to its own build tree
# (its CMakeLists' comment: PROJECT_BINARY_DIR, not CMAKE_BINARY_DIR, so it
# does not leak into a consumer's shared modules/), and exposes that
# directory to whoever links pic-mpi::pic-mpi. What it does not do is add
# `pic`'s own module directory to ITS OWN compile flags -- pic::pic here is
# the bare INTERFACE-IMPORTED alias sample_utils.cmake's my_fetch_package
# creates when nothing has published one already, and that alias carries
# no INTERFACE_INCLUDE_DIRECTORIES of its own. That is invisible in
# pic-mpi's own standalone build, where PROJECT_BINARY_DIR and
# CMAKE_BINARY_DIR coincide and `pic` writes its modules right where
# pic-mpi looks; it is not invisible here, where `pic` was already fetched
# (by pic-blas, above) into terco's shared modules/ before pic-mpi's
# CMakeLists ever runs. Without this, pic-mpi's own `use pic_types` fails
# to compile the first time TERCO_ENABLE_MPI is turned on.
if(TARGET pic-mpi AND TARGET pic)
  get_target_property(_pic_moddir pic Fortran_MODULE_DIRECTORY)
  if(_pic_moddir)
    target_include_directories(pic-mpi PRIVATE "${_pic_moddir}")
  endif()
  unset(_pic_moddir)
endif()
