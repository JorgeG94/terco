target_sources(
  ${main_lib}
  PRIVATE
          trc_scf_settings.F90       # trc_scf_options_t/result_t and the options a run resolves them to
          trc_scf_context.F90        # everything an SCF run builds once per geometry
          trc_scf_state.F90          # the SCF's per-iteration arrays, and the one enter/exit pair for them
          trc_scf_fock.F90           # the Fock build: two-electron, XC, PCM, energy, commutator, folds
          trc_scf_df.F90             # density-fitted J/K (--df): metric, Cholesky, (mn|P), B, per-iteration J/K
          trc_scf_df_response.F90    # G_fit[X] C_occ for any symmetric X: the --df reference's response potential
          trc_scf_diis.F90           # CDIIS and ADIIS, the level shift and the damped step
          trc_scf_incremental.F90    # when the Fock build is an increment and when it is full
          trc_scf_driver.F90         # the SCF: HF and Kohn-Sham, R and U, resident
          trc_scf_file.F90           # a converged SCF to a file and back, for a restart
          trc_scf_gdm.F90             # geometric direct minimisation core: L-BFGS history, preconditioner, step control
          trc_scf_gdm_step.F90        # GDM as a step mode of the SCF loop: accept-on-descent orbital rotations
          trc_scf_soscf.F90           # second-order SCF (--soscf): trust-region Newton on the orbital rotations, RHF
)

add_subdirectory(gradient)
add_subdirectory(response)
add_subdirectory(vibration)
add_subdirectory(guess)
add_subdirectory(analysis)
