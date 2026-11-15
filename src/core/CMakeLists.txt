target_sources(
  ${main_lib}
  PRIVATE trc_boys.F90       # Boys function, tabulated and asymptotic
          trc_tables.F90     # Cartesian index tables
          trc_cart.F90
          trc_ecp_quadrature.F90 # ECP radial quadrature and scaled spherical Bessel function
          trc_basis_types.F90 # trc_basis_t and trc_pairlist_t: the basis and the screened pair list
          trc_device_binding.F90 # one GPU per rank: trc_bind_device and friends
          trc_linalg.F90     # dense linear algebra on device-resident matrices
          trc_si.F90         # a quantity printed with its SI prefix (2.45 PFLOP)
          trc_step_table.F90 # shared timed-step table (RI-MP2, DF-JK): time/FLOP/B/s, NVTX
          trc_mpi_device.F90 # whether MPI gets device buffers directly (CUDA-aware)
          trc_constants.F90  # the single source of truth for math/CODATA constants and atomic masses
          trc_cart_perm.F90  # AO-index permutation, lexicographic <-> this build's order
          trc_ao_norm.F90    # AO-index normalisation, this build's <-> GAMESS's per-component convention
          trc_sph.F90        # the spherical (pure) AO space: block-sparse Cartesian <-> spherical transforms
          trc_acc_shim.F90   # see the file: this one is not optional
          trc_error.F90)
