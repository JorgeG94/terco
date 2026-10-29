"""McMurchie-Davidson core shared by the J engine, K and gradients.

`hermite` holds the Hermite index sets and the E-coefficient recursion,
`rtree` the R^(m) recursion as a pruned DAG with Palethorpe-Barca batching.
Both are stdlib only, so CMake can run the generators on a bare Python.
`reference` is the numpy/pyscf oracle and is never imported by a generator.

REFERENCES
----------
Inspiration only; the formulation and the code are terco's own.

  L. E. McMurchie, E. R. Davidson, J. Comput. Phys. 26, 218 (1978).
  E. Palethorpe, G. M. J. Barca, J. Chem. Theory Comput. 21, 9388 (2025).
  A. Asadchev, E. F. Valeev, J. Phys. Chem. A 127, 10889 (2023);
    J. Chem. Phys. 160, 244109 (2024).
"""
