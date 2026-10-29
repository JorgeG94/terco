# terco

<img src="docs/assets/powered_by_Fortran_switch_LC.svg" alt="powered by Fortran" width="320">

```
Adjective, from Spanish: Stubborn
```

Terco is a project born out of my stubborness and refusal to accept that something
cannot be done in plain old Fortran. The idea is simple: `do concurrent` is
standard Fortran and exposes parallelism nicely via a defined contract. 

Question: can we evaluate ERIs on the GPU just using do concurrent and 
have good performance? 

The question was partially answered by Melisa Alkan and collaborators
who used OpenMP target and do concurrent to parallelize the ERIs inside
GAMESS. See the [first paper](https://pubs.aip.org/aip/jcp/article/161/8/082501/3309322), 
the [second paper](https://pubs.aip.org/aip/jcp/article/164/21/214118/3393751) 
by Daniel del Angel who tested it for F functions, but limited himself to OpenMP target.

However, the algorithm used to evaluate the integrals relied on the more 
traditional per-shell-quartet approach which can quickly hit a bottleneck 
by starving the GPU a bit of work. Thus, I adopted here the batching idea
of preparing batches of shell pairs and distributing them to the device. 


Thus I think this is a novel-ish implementation which ended up being 
quite nice to work with and extend. 

The idea is simple, we use OpenACC to map the data between the host and 
the device (this wouldn't be needed in a unified memory environment), and 
rely exclusively on do concurrent to do the parallelization of the quartets. 



---

## What it does

| | |
|---|---|
| **Four-centre ERIs and Fock builds** | Head-Gordon–Pople through **d**, chunked Rys through **f** (`TERCO_LMAX=3`) and the low-L **g** classes, McMurchie–Davidson exchange for the rest of **g** and a McMurchie–Davidson J engine (`TERCO_LMAX=4`) |
| | single density, **batched over N densities in one integral pass**, non-symmetric density, and a device-resident entry point |
| | `j_scale` / `k_scale`, so a hybrid functional or a response solve needs no second routine |
| **Screening** | Schwarz bounds and optional density weighting, over a binned shell-pair container |
| **One electron** | overlap, kinetic, nuclear attraction |
| **Multipoles** | full Cartesian tensors through the **octupole** — 3 + 9 + 27 components |
| **Density fitting** | two- and three-centre Coulomb integrals |
| **Kohn-Sham DFT** | LDA, GGA and meta-GGA (tau-dependent: r2SCAN, SCAN, TPSS, the Minnesota M05-M15 family, ...) functionals, global, range-separated and double hybrids and VV10, RKS and UKS; every kernel checked against libxc ([docs/functionals.md](docs/functionals.md)) |
| **Correlation methods** | Terco runs ri-mp2 and also supports Double-Hybrid DFT energies  |
| **Solvation** | C-PCM and IEF-PCM SCF energies (pyscf's SWIG cavity, term for term), RKS and UKS, any functional, CPU and GPU |
| **Analytic gradients** | HF/UHF, LDA/GGA/hybrid/range-separated KS-DFT with VV10, and the `--df` fitted RHF/RKS gradient; validated against pyscf and finite differences; single GPU (a multi-GPU run computes the gradient on rank 0 and warns); no meta-GGA, PCM, RI-MP2 or double-hybrid gradients |
| **Interfaces** | a Fortran API, and a `bind(c)` ABI for host codes built with another compiler |

terco **ships no basis set**. A consumer passes exponents, coefficients and
centres — or hands over libcint's `atm`/`bas`/`env` directly — and terco builds
the basis object, the shell-pair container and the screening tables from that.

Terco will read the qcshema basis set format so if you have those, you can use them.

## Caveats (?) 

**The code has not been tested on AMD or Intel GPUs**. This is future work, if anyone is happy to try.

- **g functions** need `TERCO_LMAX=4` and are for energies (range-separated
  hybrids included): the non-symmetric Fock build and gradients refuse a g
  basis, and the multipole integrals stop at d. An LMAX=4 build compiles for about 50 minutes.

Spherical harmonics and range separation are supported:

- a spherical basis is computed in Cartesians and transformed, and `--cart`
  keeps a basis Cartesian;
- range-separated hybrids use an erf-attenuated `erf(ωr₁₂)/r₁₂` exchange
  kernel.

A host code is expected to check before it commits: `trc_basis_maxl` answers
before anything expensive is built.

---

## Building

The **GPU** build requires **NVHPC**. `-stdpar=gpu` is the whole mechanism —
it is what turns `do concurrent` into device kernels — so `TERCO_ENABLE_GPU=ON`
with any other compiler is refused rather than silently producing a slow CPU
library.

A **host** build (`-DTERCO_ENABLE_GPU=OFF`) works with **gfortran 15 or
newer**, ifx, or nvfortran. gfortran 15 is a floor rather than a preference:
`do concurrent` locality specifiers are Fortran 2018 and landed there. The host
build is not a fallback to run production on — it is how the library is
validated without a GPU, and it must give the same integrals, because a number
that depends on the offload model is a bug in the offload model.

In theory `do concurrent` can be parallelised over CPU threads: ifx, flang and
nvfortran will do it. This is not tested in CI at the moment. Work is underway
in gfortran to enable it.

The parallelism model exposed in terco is suited for GPU execution it might be 
crap on the CPU. Don't quote me on this (Jorge).

### CMake

```sh
cmake -S . -B build -DCMAKE_Fortran_COMPILER=nvfortran -DTERCO_GPU_ARCH=cc70
cmake --build build -j
ctest --test-dir build --output-on-failure
```

| Option | Default | |
|---|---|---|
| `TERCO_LMAX` | `2` | four-centre ceiling. Anything but 2 regenerates the kernels (needs Python 3) |
| `TERCO_GPU_ARCH` | `cc70` | as nvfortran spells it |
| `TERCO_BUILD_SHARED` | `ON` | `libterco.so`, which a host built with another compiler can link; `OFF` for a static library |
| `TERCO_REGENERATE` | `OFF` | refresh the committed kernels in-tree |
| `TERCO_FORCE_REGENERATE` | `OFF` | regenerate the build-tree kernels once, ignoring the generation stamps (a configure otherwise reruns a generator only when a generator, its inputs or its options changed) |
| `TERCO_ENABLE_MPI` | `OFF` | split the Fock build over MPI ranks, one device each; `OFF` is the same code on pic-mpi's single-rank backend |
| `TERCO_ENABLE_TESTING` | `ON` | |

Consume it with `find_package(terco)` and link `terco::terco`.

### fpm

```sh
tools/fpm_generate.sh                     # the Rys kernels, into src/generated/fpm/
fpm build --compiler mpifort --flag "-ffree-line-length-none"            # gfortran
fpm build --compiler mpifort --flag "-stdpar=gpu -acc -gpu=cc70,mem:separate"  # nvfortran
```

fpm builds the **shipped configuration and only that one**: `TRC_LMAX=2`. fpm
cannot run a code generator and most of this library is generated, so
`src/generated/` holds the committed HGP kernels, and `tools/fpm_generate.sh`
generates the Rys kernels (energy and gradient) that CMake would generate at
configure time. Rerun it after changing `scripts/`. Any other ceiling needs
CMake.

---

## Using it

### From the command line

The CMake build also makes a `terco` program (`-DTERCO_BUILD_APP=OFF` skips
it). It reads an xyz file and runs an SCF, in the spirit of xtb:

```sh
terco water.xyz                              # RHF/def2-SVP from a SAD guess
terco water.xyz -b cc-pvdz -f pbe0           # RKS/PBE0
terco oh.xyz -u 1 -f b3lyp                   # UKS doublet
terco nh4.xyz -c 1                           # cation: SADQ guess by default
terco water.xyz -b 6-31g --rimp2             # RHF + RI-MP2 (cc-pvdz-rifit)
terco water.xyz --json result.json           # also write a JSON summary
terco water.xyz -f b3lyp --pcm cpcm --eps 78.3553   # in water, C-PCM
terco --help
```

A basis name is looked up as `<name>.json` in each directory of
`$TERCO_BASIS_PATH`, then in `share/terco/basis_sets` under the install prefix,
where `cmake --install` puts the basis sets terco is tested with. The exit
status is 0 on success, 1 on an error and 2 if the SCF did not converge.

#### Over MPI ranks

Build with `-DTERCO_ENABLE_MPI=ON` (NVHPC's own MPI; see the CMake table
above), then run under `mpiexec`:

```sh
mpiexec -np 4 terco water.xyz -b 6-31g*       # one GPU per rank
```

Work is divided among GPUs in a mostly load balanced way - no dynamic
load balancing is enabled (yet). The code is able to run in a distributed
fashion, i.e. (with slurm because it is very nice to read), you do:
`srun -N 2 --ntasks=8 --ntasks-per-node=4 --ngpus-per-node=4 ./terco big.xyz`
it will run on two nodes and distribute among 8 GPUs. All methods that can 
do multi GPU can do this. I.e. just the gradients cannot.

### From Fortran

```fortran
use trc_api,  only: trc_basis_t, trc_pairlist_t, trc_1e, trc_multipoles
use trc_fock, only: trc_eri_t

type(trc_basis_t)    :: bas
type(trc_pairlist_t) :: pl
type(trc_eri_t)      :: eri

call bas%build(nsh, sh_l, sh_np, sh_e, sh_c, sh_r, nat, at_z, at_r, maxnp)
call bas%to_device()
call pl%build(bas, 1.0e-10_dp)
call pl%to_device()
call eri%build(bas, 1.0e-10_dp)          ! Schwarz bounds, bins, pair tables

call trc_1e(bas, pl, smat, tmat, vmat)   ! S, T, V in one pass
call eri%fock(bas, dmat, gmat)           ! G = J - K/2
call eri%fock_many(bas, n, dmats, gmats) ! N densities, ONE integral pass
```

The basis, the pair list and the ERI object are built **once per geometry** and
reused by every Fock build. Rebuilding them per SCF iteration throws away the
container's entire reason for existing. If used in a fragmentation code, this can 
be done better. 

`fock_many` takes `(ndens, nao, nao)` — the density index **first**. That is not
arbitrary: with it last, the N atomic updates for one `(mu,nu)` sit `nao²` apart
and the batch stops paying past N = 4.

### From another compiler

A gfortran or ifx host cannot read nvfortran's `.mod` files, so it goes through
the C ABI. `libterco.so` is the default build — nvfortran links
it, so it carries the NVHPC and OpenACC runtimes as its own `DT_NEEDED` and the
host's link line needs nothing but `-lterco`.

Compile `include/trc_c_interfaces.f90` with your own compiler and `use` it.
That file is the **authoritative declaration**; do not write your own interface
block. A copy of a contract kept in step by hand still compiles and still links
when it drifts, because a C symbol carries no signature — it fails at run time
by reading whatever was in the register.

```fortran
use trc_c_interfaces
rc = trc_basis_create_libcint(atm, natm, bas, nbas, env, nenv, 1, hbas)
rc = trc_eri_create(hbas, 1.0e-10_dp, heri)
rc = trc_fock(heri, hbas, d, g, 1.0d0, 1.0d0, 1, 0.0d0)   ! last: omega, 0 = Coulomb
```

---

## How it is put together

Most of the compiled lines come out of `scripts/`. That is the design: a
generated kernel has every index as a literal, no table loads and no branches,
and NVHPC budgets registers for one angular-momentum class rather than for the
worst class across all of them.

- **[`docs/GENERATORS.md`](docs/GENERATORS.md)** — what each generator emits,
  the conventions they share, how to extend them, and the traps that have
  actually cost time.
- **[`docs/ARCHITECTURE.md`](docs/ARCHITECTURE.md)** — the object model, the
  screening and binning design, and the density-fitting plan.

## References

Several design decisions here follow the GPU Hartree–Fock work in **[1]** — the
shell-pair binning by angular momentum and contraction depth, the pre-screen
before binning, the folded six-update digestion, and the per-order Boys
interpolation. The sources cite it as `[1]` where it is relevant. The rest is
the surrounding literature on GPU integral evaluation and on offloading GAMESS,
which is where the conventions and the comparison points come from.

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

<details>
<summary>BibTeX for all of the above</summary>

```bibtex
@article{galvez_vallejo_high-performance_2022,
  author  = {Galvez Vallejo, Jorge Luis and Barca, Giuseppe M. J. and Gordon, Mark S.},
  title   = {High-performance {GPU}-accelerated evaluation of electron repulsion integrals},
  journal = {Molecular Physics},
  volume  = {121},
  number  = {9--10},
  pages   = {e2112987},
  year    = {2022},
  issn    = {0026-8976},
  doi     = {10.1080/00268976.2022.2112987},
}

@article{del_angel_cruz_f_type_2026,
  author  = {Del Angel Cruz, Daniel and Galvez Vallejo, Jorge L. and Gordon, Mark S.},
  title   = {Electron repulsion integral evaluation over f-type functions on {GPUs} via {OpenMP} offloading},
  journal = {The Journal of Chemical Physics},
  volume  = {164},
  number  = {21},
  pages   = {214118},
  year    = {2026},
  doi     = {10.1063/5.0334077},
}

@article{alkan_liberi_2024,
  author  = {Alkan, Melisa and Pham, Buu Q. and Del Angel Cruz, Daniel and Hammond, Jeff R. and Barnes, Taylor A. and Gordon, Mark S.},
  title   = {{LibERI}---A portable and performant multi-{GPU} accelerated library for electron repulsion integrals via {OpenMP} offloading and standard language parallelism},
  journal = {The Journal of Chemical Physics},
  volume  = {161},
  number  = {8},
  pages   = {082501},
  year    = {2024},
  doi     = {10.1063/5.0215352},
}

@inproceedings{sattasathuchana_performance_portability_2024,
  author    = {Sattasathuchana, Tosaporn and Xu, Peng and Oryspayev, Dossay and Bertoni, Colleen and Roskop, Luke B. and Gordon, Mark S.},
  title     = {Performance {Portability} of {Electron} {Repulsion} {Integrals} and {Their} {Related} {Methods} across {Peta} to {Exascale} {Architectures}},
  booktitle = {SC24-W: Workshops of the International Conference for High Performance Computing, Networking, Storage and Analysis},
  publisher = {IEEE},
  pages     = {1943--1954},
  year      = {2024},
  doi       = {10.1109/SCW63240.2024.00244},
}

@article{sosonkina_runtime_2024,
  author  = {Sosonkina, Masha and Mateescu, Gabriel and Xu, Peng and Sattasathuchana, Tosaporn and Pham, Buu and Gordon, Mark S. and Leang, Sarom S.},
  title   = {Runtime performance of a {GAMESS} quantum chemistry application offloaded to {GPUs}},
  journal = {Concurrency and Computation: Practice and Experience},
  volume  = {36},
  number  = {23},
  pages   = {e8244},
  year    = {2024},
  doi     = {10.1002/cpe.8244},
}

@incollection{chapman_openmp_hackathon_2021,
  author    = {Chapman, Barbara and Pham, Buu and Yang, Charlene and Daley, Christopher and Bertoni, Colleen and Kulkarni, Dhruva and Oryspayev, Dossay and D'Azevedo, Ed and Doerfert, Johannes and Zhou, Keren and Ravikumar, Kiran and Gordon, Mark and Del Ben, Mauro and Lin, Meifeng and Alkan, Melisa and Kruse, Michael and Hernandez, Oscar and Yeung, P. K. and Lin, Paul and Xu, Peng and Pophale, Swaroop and Sattasathuchana, Tosaporn and Kale, Vivek and Huhn, William and He, Yun},
  title     = {Outcomes of {OpenMP} {Hackathon}: {OpenMP} {Application} {Experiences} with the {Offloading} {Model} ({Part} {II})},
  booktitle = {OpenMP: Enabling Massive Node-Level Parallelism},
  editor    = {McIntosh-Smith, Simon and de Supinski, Bronis R. and Klinkenberg, Jannis},
  series    = {Lecture Notes in Computer Science},
  volume    = {12870},
  publisher = {Springer International Publishing},
  pages     = {81--95},
  year      = {2021},
  isbn      = {978-3-030-85261-0},
  doi       = {10.1007/978-3-030-85262-7_6},
}

@article{sattasathuchana_gefp_2019,
  author  = {Sattasathuchana, Tosaporn and Xu, Peng and Gordon, Mark S.},
  title   = {An {Accurate} {Quantum}-{Based} {Approach} to {Explicit} {Solvent} {Effects}: {Interfacing} the {General} {Effective} {Fragment} {Potential} {Method} with \emph{Ab Initio} {Electronic} {Structure} {Theory}},
  journal = {The Journal of Physical Chemistry A},
  volume  = {123},
  number  = {39},
  pages   = {8460--8475},
  year    = {2019},
  doi     = {10.1021/acs.jpca.9b05801},
}
```

</details>

Validation and comparison targets used during development, all open source:
 [libcint], [pyscf], [gpu4pyscf], and [GAMESS libERI].

[libcint]: https://github.com/sunqm/libcint
[pyscf]: https://github.com/pyscf/pyscf
[gpu4pyscf]: https://github.com/pyscf/gpu4pyscf
[GAMESS libERI]: https://github.com/gms-bbg/gms_libERI

## Licence

MIT. See [LICENSE](LICENSE).

### Third-party

The XC kernels under `src/xc/trc_xc_*.F90` whose header says "GENERATED by
tools/exchcxx_to_fortran.py" are a line-for-line Fortran translation of
[ExchCXX](https://github.com/wavefunction91/ExchCXX)'s built-in kernels
(BSD 3-clause; Copyright (c) 2020-2024, The Regents of the University of
California, through Lawrence Berkeley National Laboratory; portions
Copyright (c) Microsoft Corporation). Each translated file's own header
carries that notice; ExchCXX's original is at
`include/exchcxx/impl/builtin/kernels/` in its own repository.

The "powered by Fortran" logo at the top of this file,
`docs/assets/powered_by_Fortran_switch_LC.svg`, is the switch logo from
[powered_by_fortran_logo](https://github.com/sebastian-mutz/powered_by_fortran_logo)
(MIT; Copyright (c) 2026 Sebastian G. Mutz), unchanged; its licence is
`docs/assets/LICENSE.powered_by_fortran_logo`.

[libxc](https://github.com/libxc/libxc) is used only as the test reference
for `check_xc` and friends, reached through [pyscf] -- never linked into
terco.
