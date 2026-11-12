# Overview

Terco is a GPU accelerated quantum chemistry library written in standard Fortran.

The goal of terco is to showcase modern Fortran and GPU performance without writing
a single line of CUDA and relying on OpenACC/OpenMP for data transfers. The second goal
is to motivate compiler developers across the vendors to solidify `do concurrent` support
across hardware, since in theory, it should be portable.

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
- **Correlation**: RI-MP2, and double-hybrid DFT energies; closed-shell
  RI-CCSD and RI-CCSD(T) energies (RHF reference, no gradients), on one GPU or
  over several with one MPI rank per GPU. See {doc}`theory` for the method and
  its limits, and {doc}`cli` for the options.
- **Solvation**: C-PCM and IEF-PCM SCF energies, RKS and UKS, any
  functional, CPU and GPU.
- **Analytic gradients**: HF/UHF and KS-DFT (LDA/GGA/meta-GGA/hybrid/
  range-separated) with VV10, the density-fitted RHF/RKS gradient, RHF +
  RI-MP2, and double hybrids; four-centre or density-fitted reference;
  frozen core; Cartesian or spherical; restricted; single GPU.
- **Point-group symmetry** (opt-in, `--symmetry`): Abelian groups (D2h and
  subgroups) detected and symmetrised, an SCF blocked by irrep with orbital
  labels, and a skeleton Fock build and gradient from orbit representatives of
  the integral quartets; see {doc}`theory`.
- **Interfaces**: a Fortran API (`trc_api`), a `bind(c)` ABI for callers
  built with another compiler, and the `terco` command-line program.
- **Geometry optimization**: Through LIBDLFIND.

The library takes its basis from the caller: exponents, coefficients and
centres, or libcint's `atm`/`bas`/`env` layout directly, and builds the
basis object, the shell-pair container and the screening tables from that.
BSE-format JSON basis files and the QCSchema basis format are both read
directly, and the `terco` program ships a library of Basis Set Exchange sets
(see {doc}`install`).

## Notes on MPI

The Fock build is naively distributed, i.e. I did no galaxy brain techniques to reduce 
communication, optimize for size or anything, it is a good ol' every rank has a copy of
the Fock matrix, all reduce at the end of each iteraiton, and the diagonalization 
runs on a single GPU because ain't no-one is doing multi-GPU diagonalization because
it does not scale. Open issue: implement diag free approaches, I welcome pull requests.

## Caveats

- The code has not been validated on AMD or Intel GPUs.
- **g functions** need `TERCO_LMAX=4` and are for energies only (range-
  separated hybrids included): the non-symmetric Fock build and gradients
  refuse a g basis, and the multipole integrals stop at d.
- **RI-CCSD(T)** is closed-shell and RHF only, energies only, with no point-group
  symmetry (a named group is refused), no effective core potential and an aufbau
  reference; its size is bounded by the memory of one GPU ({doc}`theory`).
- `trc_basis_maxl` lets a host code check the highest angular momentum in a
  basis before it commits to building anything expensive with it.

See the top-level [README](https://github.com/JorgeG94/terco#readme) for the
full capability table, the references this work builds on, and the licence.
