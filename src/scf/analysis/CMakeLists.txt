target_sources(
  ${main_lib}
  PRIVATE
          trc_hirshfeld_data.F90  # the Hirshfeld proatom radial tables, loaded from a data file
          trc_hirshfeld.F90    # Hirshfeld atomic charges: stockholder partition of the molecular density
          trc_fukui.F90        # Hirshfeld-condensed Fukui indices: f+/f-/f0 from three states' charges
          trc_chelpg.F90       # CHELPG atomic charges: point charges fitted to the electrostatic potential on a grid
)
