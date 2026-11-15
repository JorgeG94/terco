target_sources(
  ${main_lib}
  PRIVATE
          trc_eri_data.F90    # the containers the four-centre launchers share
          trc_batch.F90       # primitive-pair batching
          trc_hgp.F90         # Head-Gordon-Pople driver
          trc_screen.F90      # Schwarz bounds
          trc_parity.F90      # parity skipping: dead axis, shell mask, axis permutation
          trc_bins.F90        # the binned shell-pair container
          trc_bin_tiles.F90   # exchange tiles: density pruning of the launch plan
          trc_binkernel.F90   # the shared quartet kernel and the digestion
          trc_grad.F90        # two-electron gradient launch plan, Rys only
          trc_integrals.F90   # the 1e, multipole, PCM and DF integral drivers
          trc_ecp_angular.F90  # ECP angular factors (type 1 and type 2)
          trc_ecp_radial.F90   # ECP radial factors (type 1 and type 2)
          trc_ecp.F90           # the scalar ECP one-electron matrix
          trc_mdj.F90        # the McMurchie-Davidson J engine (no quartets)
          trc_mdg.F90        # MD exchange for the g classes above TERCO_RYS_G_LTOT
          trc_eri.F90        # the four-centre Fock builds
)
