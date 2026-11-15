# Lebedev spheres, Treutler-Ahlrichs radial mesh, pruning, and partitions.
target_sources(
  ${main_lib}
  PRIVATE
          trc_lebedev_data.F90
          trc_lebedev.F90
          trc_dft_radial_data.F90
          trc_dft_radial.F90
          trc_dft_partition.F90
          trc_dft_prune.F90
          trc_dft_grid.F90
)
