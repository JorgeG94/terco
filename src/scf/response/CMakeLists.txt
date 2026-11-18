target_sources(
  ${main_lib}
  PRIVATE
          trc_cphf.F90                # coupled-perturbed HF: the orbital-response solve and the static polarizability
          trc_cpks.F90                # coupled-perturbed KS: trc_cphf plus the XC kernel response, RKS LDA/GGA
          trc_tddft.F90               # TDA/CIS and full-RPA singlet excitation energies, block Davidson over fock_nosym_many + fxc_apply_many
          trc_tddft_grad.F90          # TDA/CIS excited-state nuclear gradients for several roots: Z-vector, one gradient_pairs sweep, XC response
          trc_tddft_rpa_grad.F90      # the same for the full RPA (TDHF / Casida TDDFT): X + Y and X - Y in the Z-vector and the derivative pairs
          trc_tddft_run.F90           # the entry the C ABI and CLI call: energies and one root's gradient from a converged SCF, grid response on by default
)
