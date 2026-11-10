# Citing terco

If you use terco, please cite **[1]** below, the GPU Hartree–Fock paper several of terco's
design decisions follow (the repository's `CITATION.cff` says the same). When a feature
rests on a method listed further down, please cite that too. Each of these is cited in
terco's source where it is used, and all of them are in
{download}`citations.bib <citations.bib>` (generated from Crossref records).

## terco's background

The papers listed in the README: the GPU Hartree–Fock work in **[1]** (shell-pair binning
by angular momentum and contraction depth, the pre-screen before binning, the folded
six-update digestion, the per-order Boys interpolation), and the surrounding work on GPU
integrals and on offloading GAMESS.

> **[1]** J. L. Galvez Vallejo, G. M. J. Barca and M. S. Gordon,
> *High-performance GPU-accelerated evaluation of electron repulsion integrals*,
> **Molecular Physics** **121** (2022), e2112987.
> [doi:10.1080/00268976.2022.2112987](https://doi.org/10.1080/00268976.2022.2112987)

> **[2]** D. Del Angel Cruz, J. L. Galvez Vallejo and M. S. Gordon,
> *Electron repulsion integral evaluation over f-type functions on GPUs via
> OpenMP offloading*,
> **The Journal of Chemical Physics** **164** (2026), 214118.
> [doi:10.1063/5.0334077](https://doi.org/10.1063/5.0334077)

> **[3]** M. Alkan, B. Q. Pham, D. Del Angel Cruz, J. R. Hammond, T. A. Barnes
> and M. S. Gordon,
> *LibERI — A portable and performant multi-GPU accelerated library for electron
> repulsion integrals via OpenMP offloading and standard language parallelism*,
> **The Journal of Chemical Physics** **161** (2024), 082501.
> [doi:10.1063/5.0215352](https://doi.org/10.1063/5.0215352)

> **[4]** T. Sattasathuchana, P. Xu, D. Oryspayev, C. Bertoni, L. B. Roskop and
> M. S. Gordon,
> *Performance portability of electron repulsion integrals and their related
> methods across peta to exascale architectures*,
> in **SC24-W: Workshops of the International Conference for High Performance
> Computing, Networking, Storage and Analysis**, IEEE (2024), 1943–1954.
> [doi:10.1109/SCW63240.2024.00244](https://doi.org/10.1109/SCW63240.2024.00244)

> **[5]** M. Sosonkina, G. Mateescu, P. Xu, T. Sattasathuchana, B. Pham,
> M. S. Gordon and S. S. Leang,
> *Runtime performance of a GAMESS quantum chemistry application offloaded to
> GPUs*,
> **Concurrency and Computation: Practice and Experience** **36** (2024), e8244.
> [doi:10.1002/cpe.8244](https://doi.org/10.1002/cpe.8244)

> **[6]** B. Chapman, B. Pham, C. Yang, C. Daley, C. Bertoni, D. Kulkarni,
> D. Oryspayev, E. D'Azevedo, J. Doerfert, K. Zhou, K. Ravikumar, M. Gordon,
> M. Del Ben, M. Lin, M. Alkan, M. Kruse, O. Hernandez, P. K. Yeung, P. Lin,
> P. Xu, S. Pophale, T. Sattasathuchana, V. Kale, W. Huhn and Y. He,
> *Outcomes of OpenMP hackathon: OpenMP application experiences with the
> offloading model (Part II)*,
> in **OpenMP: Enabling Massive Node-Level Parallelism**, edited by
> S. McIntosh-Smith, B. R. de Supinski and J. Klinkenberg, Lecture Notes in
> Computer Science **12870**, Springer (2021), 81–95.
> [doi:10.1007/978-3-030-85262-7_6](https://doi.org/10.1007/978-3-030-85262-7_6)

> **[7]** T. Sattasathuchana, P. Xu and M. S. Gordon,
> *An accurate quantum-based approach to explicit solvent effects: interfacing
> the general effective fragment potential method with ab initio electronic
> structure theory*,
> **The Journal of Physical Chemistry A** **123** (2019), 8460–8475.
> [doi:10.1021/acs.jpca.9b05801](https://doi.org/10.1021/acs.jpca.9b05801)

## Integral algorithms

The recurrences and quadratures terco's ERI kernels are built on.

- L. E. McMurchie and E. R. Davidson, *One- and two-electron integrals over cartesian gaussian functions*, Journal of Computational Physics **26**, 218–231 (1978). [doi:10.1016/0021-9991(78)90092-x](https://doi.org/10.1016/0021-9991(78)90092-x)
- S. Obara and A. Saika, *Efficient recursive computation of molecular integrals over Cartesian Gaussian functions*, The Journal of Chemical Physics **84**, 3963–3974 (1986). [doi:10.1063/1.450106](https://doi.org/10.1063/1.450106)
- M. Head-Gordon and J. A. Pople, *A method for two-electron Gaussian integral and integral derivative evaluation using recurrence relations*, The Journal of Chemical Physics **89**, 5777–5786 (1988). [doi:10.1063/1.455553](https://doi.org/10.1063/1.455553)
- M. Dupuis, J. Rys and H. F. King, *Evaluation of molecular integrals over Gaussian basis functions*, The Journal of Chemical Physics **65**, 111–116 (1976). [doi:10.1063/1.432807](https://doi.org/10.1063/1.432807)

## GPU integral and Fock-build work terco learned from

Inspiration: the formulation and the code in terco are its own, and each source file names what it drew on.

- I. S. Ufimtsev and T. J. Martínez, *Quantum Chemistry on Graphical Processing Units. 1. Strategies for Two-Electron Integral Evaluation*, Journal of Chemical Theory and Computation **4**, 222–231 (2008). [doi:10.1021/ct700268q](https://doi.org/10.1021/ct700268q)
- I. S. Ufimtsev and T. J. Martínez, *Quantum Chemistry on Graphical Processing Units. 2. Direct Self-Consistent-Field Implementation*, Journal of Chemical Theory and Computation **5**, 1004–1015 (2009). [doi:10.1021/ct800526s](https://doi.org/10.1021/ct800526s)
- Y. Wang, D. Hait, K. G. Johnson, O. J. Fajen, J. H. Zhang, R. D. Guerrero and T. J. Martínez, *Extending GPU-accelerated Gaussian integrals in the TeraChem software package to f type orbitals: Implementation and applications*, The Journal of Chemical Physics **161**, 174118 (2024). [doi:10.1063/5.0233523](https://doi.org/10.1063/5.0233523)
- A. Asadchev and E. F. Valeev, *High-Performance Evaluation of High Angular Momentum 4-Center Gaussian Integrals on Modern Accelerated Processors*, The Journal of Physical Chemistry A **127**, 10889–10895 (2023). [doi:10.1021/acs.jpca.3c04574](https://doi.org/10.1021/acs.jpca.3c04574)
- A. Asadchev and E. F. Valeev, *3-center and 4-center 2-particle Gaussian AO integrals on modern accelerated processors*, The Journal of Chemical Physics **160**, 244109 (2024). [doi:10.1063/5.0217001](https://doi.org/10.1063/5.0217001)
- R. Li, Q. Sun, X. Zhang and G. K.-L. Chan, *Introducing GPU Acceleration into the Python-Based Simulations of Chemistry Framework*, The Journal of Physical Chemistry A **129**, 1459–1468 (2025). [doi:10.1021/acs.jpca.4c05876](https://doi.org/10.1021/acs.jpca.4c05876)
- E. Palethorpe and G. M. J. Barca, *High-Performance, High-Angular-Momentum J Engine on Graphics Processing Units*, Journal of Chemical Theory and Computation **21**, 9388–9403 (2025). [doi:10.1021/acs.jctc.5c00775](https://doi.org/10.1021/acs.jctc.5c00775)
- E. Palethorpe, R. Stocks and G. M. J. Barca, *High-Performance, High-Angular-Momentum K Engine on Graphics Processing Units*, ChemRxiv (preprint) (2026). [doi:10.26434/chemrxiv.15004739/v2](https://doi.org/10.26434/chemrxiv.15004739/v2)

## Density fitting, RI-HF/RI-MP2 and DFT on GPUs

- R. Stocks, E. Palethorpe and G. M. J. Barca, *Multi-GPU RI-HF Energies and Analytic Gradients─Toward High-Throughput Ab Initio Molecular Dynamics*, Journal of Chemical Theory and Computation **20**, 7503–7515 (2024). [doi:10.1021/acs.jctc.4c00877](https://doi.org/10.1021/acs.jctc.4c00877)
- J. J. Hayes and M. S. Gordon, *Speeding Up Hartree–Fock in JuliaChem with Density Fitting*, Journal of Chemical Theory and Computation **22**, 3897–3915 (2026). [doi:10.1021/acs.jctc.5c02025](https://doi.org/10.1021/acs.jctc.5c02025)
- R. Stocks, E. Palethorpe and G. M. J. Barca, *High-Performance Multi-GPU Analytic RI-MP2 Energy Gradients*, Journal of Chemical Theory and Computation **20**, 2505–2519 (2024). [doi:10.1021/acs.jctc.3c01424](https://doi.org/10.1021/acs.jctc.3c01424)
- R. Stocks and G. M. J. Barca, *Efficient Algorithms for GPU Accelerated Evaluation of the DFT Exchange-Correlation Functional*, Journal of Chemical Theory and Computation **21**, 10263–10280 (2025). [doi:10.1021/acs.jctc.5c01229](https://doi.org/10.1021/acs.jctc.5c01229)

## SCF convergence and guesses

- P. Pulay, *Convergence acceleration of iterative sequences. The case of SCF iteration*, Chemical Physics Letters **73**, 393–398 (1980). [doi:10.1016/0009-2614(80)80396-4](https://doi.org/10.1016/0009-2614(80)80396-4)
- P. Pulay, *Improved SCF convergence acceleration*, Journal of Computational Chemistry **3**, 556–560 (1982). [doi:10.1002/jcc.540030413](https://doi.org/10.1002/jcc.540030413)
- X. Hu and W. Yang, *Accelerating self-consistent field convergence with the augmented Roothaan–Hall energy function*, The Journal of Chemical Physics **132**, 054109 (2010). [doi:10.1063/1.3304922](https://doi.org/10.1063/1.3304922)
- S. Lehtola, *Assessment of Initial Guesses for Self-Consistent Field Calculations. Superposition of Atomic Potentials: Simple yet Efficient*, Journal of Chemical Theory and Computation **15**, 1593–1604 (2019). [doi:10.1021/acs.jctc.8b01089](https://doi.org/10.1021/acs.jctc.8b01089)

## DFT grids and functionals

Every functional is checked against libxc, which lists each functional's own paper; see also {doc}`functionals`.

- A. D. Becke, *A multicenter numerical integration scheme for polyatomic molecules*, The Journal of Chemical Physics **88**, 2547–2553 (1988). [doi:10.1063/1.454033](https://doi.org/10.1063/1.454033)
- O. Treutler and R. Ahlrichs, *Efficient molecular numerical integration schemes*, The Journal of Chemical Physics **102**, 346–354 (1995). [doi:10.1063/1.469408](https://doi.org/10.1063/1.469408)
- R. E. Stratmann, G. E. Scuseria and M. J. Frisch, *Achieving linear scaling in exchange-correlation density functional quadratures*, Chemical Physics Letters **257**, 213–223 (1996). [doi:10.1016/0009-2614(96)00600-8](https://doi.org/10.1016/0009-2614(96)00600-8)
- V. I. Lebedev and D. N. Laikov, *A quadrature formula for the sphere of the 131st algebraic order of accuracy*, Doklady Mathematics **59**, 477–481 (1999).
- S. Lehtola, C. Steigemann, M. J. T. Oliveira and M. A. L. Marques, *Recent developments in libxc — A comprehensive library of functionals for density functional theory*, SoftwareX **7**, 1–5 (2018). [doi:10.1016/j.softx.2017.11.002](https://doi.org/10.1016/j.softx.2017.11.002)
- O. A. Vydrov and T. Van Voorhis, *Nonlocal van der Waals density functional: The simpler the better*, The Journal of Chemical Physics **133**, 244103 (2010). [doi:10.1063/1.3521275](https://doi.org/10.1063/1.3521275)
- N. Mardirossian and M. Head-Gordon, *ωB97M-V: A combinatorially optimized, range-separated hybrid, meta-GGA density functional with VV10 nonlocal correlation*, The Journal of Chemical Physics **144**, 214110 (2016). [doi:10.1063/1.4952647](https://doi.org/10.1063/1.4952647)
- S. Grimme, *Semiempirical hybrid density functional with perturbative second-order correlation*, The Journal of Chemical Physics **124**, 034108 (2006). [doi:10.1063/1.2148954](https://doi.org/10.1063/1.2148954)
- A. Karton, A. Tarnopolsky, J.-F. Lamère, G. C. Schatz and J. M. L. Martin, *Highly Accurate First-Principles Benchmark Data Sets for the Parametrization and Validation of Density Functional and Other Approximate Methods. Derivation of a Robust, Generally Applicable, Double-Hybrid Functional for Thermochemistry and Thermochemical Kinetics*, The Journal of Physical Chemistry A **112**, 12868–12886 (2008). [doi:10.1021/jp801805p](https://doi.org/10.1021/jp801805p)
- T. Schwabe and S. Grimme, *Towards chemical accuracy for the thermodynamics of large molecules: new hybrid density functionals including non-local correlation effects*, Physical Chemistry Chemical Physics **8**, 4398 (2006). [doi:10.1039/b608478h](https://doi.org/10.1039/b608478h)

## Response, charges and solvation

- J. A. Pople, R. Krishnan, H. B. Schlegel and J. S. Binkley, *Derivative studies in Hartree–Fock and Møller–Plesset theories*, International Journal of Quantum Chemistry **16**, 225–241 (1979). [doi:10.1002/qua.560160825](https://doi.org/10.1002/qua.560160825)
- F. L. Hirshfeld, *Bonded-atom fragments for describing molecular charge densities*, Theoretica Chimica Acta **44**, 129–138 (1977). [doi:10.1007/bf00549096](https://doi.org/10.1007/bf00549096)
- C. M. Breneman and K. B. Wiberg, *Determining atom‐centered monopoles from molecular electrostatic potentials. The need for high sampling density in formamide conformational analysis*, Journal of Computational Chemistry **11**, 361–373 (1990). [doi:10.1002/jcc.540110311](https://doi.org/10.1002/jcc.540110311)
- A. Bondi, *van der Waals Volumes and Radii*, The Journal of Physical Chemistry **68**, 441–451 (1964). [doi:10.1021/j100785a001](https://doi.org/10.1021/j100785a001)
- R. S. Rowland and R. Taylor, *Intermolecular Nonbonded Contact Distances in Organic Crystal Structures: Comparison with Distances Expected from van der Waals Radii*, The Journal of Physical Chemistry **100**, 7384–7391 (1996). [doi:10.1021/jp953141+](https://doi.org/10.1021/jp953141+)
- M. Mantina, A. C. Chamberlin, R. Valero, C. J. Cramer and D. G. Truhlar, *Consistent van der Waals Radii for the Whole Main Group*, The Journal of Physical Chemistry A **113**, 5806–5812 (2009). [doi:10.1021/jp8111556](https://doi.org/10.1021/jp8111556)
- A. W. Lange and J. M. Herbert, *A smooth, nonsingular, and faithful discretization scheme for polarizable continuum models: The switching/Gaussian approach*, The Journal of Chemical Physics **133**, 244111 (2010). [doi:10.1063/1.3511297](https://doi.org/10.1063/1.3511297)
- B. A. Gregersen and D. M. York, *High-order discretization schemes for biochemical applications of boundary element solvation and variational electrostatic projection methods*, The Journal of Chemical Physics **122**, 194110 (2005). [doi:10.1063/1.1899146](https://doi.org/10.1063/1.1899146)

## Geometry optimisation

- J. Kästner, J. M. Carr, T. W. Keal, W. Thiel, A. Wander and P. Sherwood, *DL-FIND: An Open-Source Geometry Optimizer for Atomistic Simulations*, The Journal of Physical Chemistry A **113**, 11856–11865 (2009). [doi:10.1021/jp9028968](https://doi.org/10.1021/jp9028968)

## Basis sets and reference software

terco's numbers are checked against these; the basis sets it ships come from the Basis Set Exchange.

- B. P. Pritchard, D. Altarawy, B. Didier, T. D. Gibson and T. L. Windus, *New Basis Set Exchange: An Open, Up-to-Date Resource for the Molecular Sciences Community*, Journal of Chemical Information and Modeling **59**, 4814–4820 (2019). [doi:10.1021/acs.jcim.9b00725](https://doi.org/10.1021/acs.jcim.9b00725)
- Q. Sun, *Libcint: An efficient general integral library for Gaussian basis functions*, Journal of Computational Chemistry **36**, 1664–1671 (2015). [doi:10.1002/jcc.23981](https://doi.org/10.1002/jcc.23981)
- Q. Sun, X. Zhang, S. Banerjee, P. Bao, M. Barbry, N. S. Blunt et al., *Recent developments in the PySCF program package*, The Journal of Chemical Physics **153**, 024109 (2020). [doi:10.1063/5.0006074](https://doi.org/10.1063/5.0006074)

## Software terco is built with

- [pic](https://github.com/JorgeG94/pic), [pic-blas](https://github.com/JorgeG94/pic-blas) and [pic-mpi](https://github.com/JorgeG94/pic-mpi): kinds, logging, timers, BLAS/LAPACK and MPI wrappers.
- [json-fortran](https://github.com/jacobwilliams/json-fortran) and [TOML Fortran](https://github.com/toml-f/toml-f): basis files, JSON output and the TOML input.
- [test-drive](https://github.com/fortran-lang/test-drive) for the unit tests; [FORD](https://github.com/Fortran-FOSS-Programmers/ford) and fortran-lang's [sphinx-fortran-domain](https://github.com/fortran-lang/sphinx-fortran-domain) for these docs.
- libfint and [PySCF](https://github.com/pyscf/pyscf): the ECP integrals are a Fortran translation of libfint's, itself a translation of PySCF's (Apache-2.0; see `NOTICE`).
- [metalquicha](https://github.com/JorgeG94/metalquicha): several terco paths (PCM, the RI-MP2 gradient, guesses) were ported from it, and it supplies reference values.
- [GPU4PySCF](https://github.com/pyscf/gpu4pyscf): the main performance comparison (its integral paper is Li et al. above).
