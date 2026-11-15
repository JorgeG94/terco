target_sources(
  ${main_lib}
  PRIVATE
          trc_rccsd_steps.F90         # RI-CCSD(T): the step table's rows and the flop helpers every stage shares
          trc_rccsd_ints.F90          # RI-CCSD(T): the fitted MO integral blocks (oooo, ooov, oovv, ovov, ovvv, B)
          trc_rccsd_inter.F90         # RI-CCSD(T): tau, Foo/Fvv/Fov/Loo/Lvv, the energy functional, MP2 amplitudes
          trc_rccsd_w.F90             # RI-CCSD(T): the two-particle intermediates Woooo, Wvoov, Wvovo
          trc_rccsd_ladder.F90        # RI-CCSD(T): the particle-particle ladder (assembled and direct forms) and its t1 dressing
          trc_rccsd_update.F90        # RI-CCSD(T): the T1/T2 amplitude update (Hirata Eqs. 35-36) from the modules above
          trc_rccsd_diis.F90          # RI-CCSD(T): amplitude DIIS and the step measure, on the device
          trc_rccsd_mpi.F90             # RI-CCSD(T): the collective calls of the multi-rank run
          trc_rccsd_xfer.F90           # RI-CCSD(T): amplitudes and partial ladders between ranks, packed, host or device buffers
          trc_rccsd_share.F90           # RI-CCSD(T): the d ranges of the ladder over the ranks of a run
          trc_rccsd_dist.F90            # RI-CCSD(T): the top of an iteration over ranks, and a follower's loop
          trc_rccsd_driver.F90        # RI-CCSD(T): the iteration driver, convergence, result and handle
          trc_rccsd_triples_dev.F90   # RI-CCSD(T): the (T) correction on the device, tiled gemms and one energy kernel
          trc_rccsd_triples_ijk.F90   # RI-CCSD(T): the (T) correction on the device, driven by occupied triples
          trc_rccsd_triples.F90       # RI-CCSD(T): the (T) correction, host reference (Rendell-Lee-Komornicki)
)
