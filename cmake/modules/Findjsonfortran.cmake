#
# json-fortran: the JSON reader behind trc_basis_json, which builds a basis
# from the MolSSI BSE files metalquicha ships. Fetched the way pic is.
#
set(_lib "jsonfortran")
set(_url "https://github.com/jacobwilliams/json-fortran.git")
include("${CMAKE_CURRENT_LIST_DIR}/sample_utils.cmake")
set(_rev "9.3.1")
set(SKIP_DOC_GEN ON CACHE BOOL "" FORCE)
set(JSONFORTRAN_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(JSONFORTRAN_STATIC_LIBRARY_ONLY ON CACHE BOOL "" FORCE)
my_fetch_package("${_lib}" "${_url}" "${_rev}")
# json-fortran writes its module files at the root of its binary directory
# and exports no include directory a consumer could use; say where they are.
target_include_directories(jsonfortran INTERFACE ${jsonfortran_BINARY_DIR}
                                                 ${jsonfortran_BINARY_DIR}/include)
unset(_lib)
unset(_url)
unset(_rev)
