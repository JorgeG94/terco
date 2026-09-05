#
# pic-mpi: one API over mpi_f08, the legacy MPI module and a single-rank
# serial backend, so terco depends on it unconditionally and TERCO_ENABLE_MPI
# only chooses the backend (through PIC_ENABLE_MPI). No #ifdef in terco.
#
set(_lib "pic-mpi")
set(_url "https://github.com/JorgeG94/pic-mpi/")
include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")
set(_rev "v0.6.0")
my_fetch_package("${_lib}" "${_url}" "${_rev}")
unset(_lib)
unset(_url)
unset(_rev)
