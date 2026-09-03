#
# pic-blas: explicit-interface wrappers over BLAS and LAPACK, fetched the
# way pic is. The SCF driver's host-side linear algebra goes through it so
# nothing in terco is called through an implicit interface.
#
set(_lib "pic-blas")
set(_url "https://github.com/JorgeG94/pic-blas/")
include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")
set(_rev "v0.3.1")
my_fetch_package("${_lib}" "${_url}" "${_rev}")
unset(_lib)
unset(_url)
unset(_rev)
