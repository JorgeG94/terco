target_sources(
  ${main_lib}
  PRIVATE
          trc_cphf.F90                # coupled-perturbed HF: the orbital-response solve and the static polarizability
          trc_cpks.F90                # coupled-perturbed KS: trc_cphf plus the XC kernel response, RKS LDA/GGA
          trc_tddft.F90               # TDA/CIS singlet excitation energies, block Davidson over fock_split_many + fxc_apply_many
          trc_tddft_grad.F90          # TDA/CIS excited-state nuclear gradients for several roots: Z-vector, one gradient_pairs sweep, XC response
)
