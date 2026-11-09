# Overview

terco is a GPU electron-repulsion-integral (ERI) and SCF library written in
standard Fortran. The compute kernels are `do concurrent`, compiled by
nvfortran with `-stdpar=gpu -acc -gpu=mem:separate`; OpenACC is used only to
move data between host and device, never to launch compute. The same source
also builds and validates on the CPU with gfortran or ifx -- `do concurrent`
is then an ordinary serial loop and every OpenACC directive is an inert
comment, so the host build is how the library's integrals are checked
without a GPU.

## What v0.1.0 covers

- **Four-centre ERIs and Fock builds**: Head-Gordon-Pople through **d**,
  chunked Rys quadrature through **f** (`TERCO_LMAX=3`) and the low-angular-
  momentum **g** classes, and McMurchie-Davidson for the rest of **g**
  (`TERCO_LMAX=4`) and for the density-fitting J engine. Single density,
  batched over N densities in one integral pass, non-symmetric density, and
  a device-resident entry point; `j_scale`/`k_scale` so a hybrid functional
  or a response solve needs no second routine.
- **Screening**: Schwarz bounds, with optional density weighting, over a
  binned shell-pair container.
- **One-electron integrals**: overlap, kinetic, nuclear attraction.
- **Multipoles**: full Cartesian tensors through the octupole.
- **Density fitting**: two- and three-centre Coulomb integrals.
- **Kohn-Sham DFT**: LDA, GGA and meta-GGA (including tau-dependent
  functionals), global and range-separated hybrids, double hybrids, and
  VV10; restricted and unrestricted. Every kernel is checked against libxc
  -- see {doc}`functionals` for the full, generated list.
- **Correlation**: RI-MP2, and double-hybrid DFT energies.
- **Solvation**: C-PCM and IEF-PCM SCF energies, RKS and UKS, any
  functional, CPU and GPU.
- **Analytic gradients**: HF/UHF and KS-DFT (LDA/GGA/meta-GGA/hybrid/
  range-separated) with VV10, the density-fitted RHF/RKS gradient, RHF +
  RI-MP2, and double hybrids; four-centre or density-fitted reference;
  frozen core; Cartesian or spherical; restricted; single GPU.
- **Interfaces**: a Fortran API (`trc_api`), a `bind(c)` ABI for callers
  built with another compiler, and the `terco` command-line program.

terco **ships no basis set**. A caller passes exponents, coefficients and
centres -- or hands over libcint's `atm`/`bas`/`env` layout directly -- and
terco builds the basis object, the shell-pair container and the screening
tables from that. BSE-format JSON basis files and the QCSchema basis format
are both read directly.

## Caveats

- The code has not been validated on AMD or Intel GPUs.
- **g functions** need `TERCO_LMAX=4` and are for energies only (range-
  separated hybrids included): the non-symmetric Fock build and gradients
  refuse a g basis, and the multipole integrals stop at d.
- `trc_basis_maxl` lets a host code check the highest angular momentum in a
  basis before it commits to building anything expensive with it.

See the top-level [README](https://github.com/JorgeG94/terco#readme) for the
full capability table, the references this work builds on, and the licence.
