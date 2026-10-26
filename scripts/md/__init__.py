"""McMurchie-Davidson core shared by the J engine, K and gradients.

`hermite` holds the Hermite index sets and the E-coefficient recursion,
`rtree` the R^(m) recursion as a pruned DAG with Palethorpe-Barca batching.
Both are stdlib only, so CMake can run the generators on a bare Python.
`reference` is the numpy/pyscf oracle and is never imported by a generator.
"""
