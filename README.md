# terco

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
| **Four-centre ERIs and Fock builds** | Head-Gordon–Pople, through **d** |
| | single density, **batched over N densities in one integral pass**, non-symmetric density, and a device-resident entry point |
| | `j_scale` / `k_scale`, so a hybrid functional or a response solve needs no second routine |
| **Screening** | Schwarz bounds and optional density weighting, over a binned shell-pair container |
| **One electron** | overlap, kinetic, nuclear attraction |
| **Multipoles** | full Cartesian tensors through the **octupole** — 3 + 9 + 27 components |
| **Density fitting** | two- and three-centre Coulomb integrals |
| **Interfaces** | a Fortran API, and a `bind(c)` ABI for host codes built with another compiler |

terco **ships no basis set**. A consumer passes exponents, coefficients and
centres — or hands over libcint's `atm`/`bas`/`env` directly — and terco builds
the basis object, the shell-pair container and the screening tables from that.

## What it does not do

Stated plainly, because a library that fails quietly is worse than one that
refuses:

- **Angular momentum above d, on four centres.** The four-centre kernels are
  generated one per `(la,lb,lc,ld)` class, so they grow as `L⁴`: 16 classes at
  p, 81 at d, and `(dd|dd)` alone is twenty thousand lines. f would need a
  larger Boys table as well — the current one caps `la+lb+lc+ld` at 8 and
  `(ff|ff)` needs 12.
- **Spherical harmonics.** Cartesian only. A caller passing a spherical basis
  is refused rather than silently reinterpreted, because reinterpreting it
  gives wrong numbers instead of an error.
- **Range separation.** No erf-attenuated kernel.
- **More than one GPU for the XC integration.** The Fock build splits over
  MPI ranks with one device each (`TERCO_ENABLE_MPI`, `trc_bind_device`, a
  communicator to `eri%build` or `trc_scf_run`); the exchange-correlation
  integration still runs whole on every rank.

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

### CMake

```sh
cmake -S . -B build \
  -DCMAKE_Fortran_COMPILER=nvfortran \
  -DTERCO_GPU_ARCH=cc70
cmake --build build -j
ctest --test-dir build --output-on-failure
```

| Option | Default | |
|---|---|---|
| `TERCO_LMAX` | `2` | four-centre ceiling. Anything but 2 regenerates the kernels (needs Python 3) |
| `TERCO_GPU_ARCH` | `cc70` | as nvfortran spells it |
| `TERCO_BUILD_SHARED` | `ON` | `libterco.so`, which a host built with another compiler can link; `OFF` for a static library |
| `TERCO_REGENERATE` | `OFF` | refresh the committed kernels in-tree |
| `TERCO_ENABLE_MPI` | `OFF` | split the Fock build over MPI ranks, one device each; `OFF` is the same code on pic-mpi's single-rank backend |
| `TERCO_ENABLE_TESTING` | `ON` | |

Consume it with `find_package(terco)` and link `terco::terco`.

### fpm

```sh
fpm build --compiler nvfortran \
  --flag "-stdpar=gpu -acc -gpu=cc70,mem:separate -Mpreprocess \
          -DTRC_LMAX=2 -DTRC_FOCK6 -DTRC_PERCLASS"
```

fpm builds the **shipped configuration and only that one**. fpm cannot run a
code generator and most of this library is generated, so `src/generated/` holds
one committed set — `TRC_LMAX=2` — and that is what fpm compiles. Any other
ceiling needs CMake.

---

## Using it

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
container's entire reason for existing.

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
rc = trc_fock(heri, hbas, d, g, 1.0d0, 1.0d0, 1)
```

This hasn't been tested thoroughly so be wary! 

---

## Validation

Nothing here asserts a number terco produced and then froze.

`check_selftest` needs **no external program and no recorded data** — closed
forms and identities only, so it runs anywhere:

| | |
|---|---|
| Overlap against the closed form for `(s\|s)` | **3.3e-16** |
| Overlap, kinetic, nuclear attraction and the Fock matrix under a **translation of every centre** | **1.7e-15** to **1.8e-14** |
| Dipole **origin shift** against the overlap, `mu(O+a) = mu(O) - a S` | **4.4e-16** |
| `G` symmetric for a symmetric density | **exact** |
| Batched build against the same densities singly | **5.6e-17** |

Translational invariance is the one that earns its place: absolute positions
do appear, in the Gaussian product centre, so it is not automatic, and it
exercises the pair list, the screening, the bins, the kernels and the
digestion in one number.

The rest compare against an independent implementation, which is the only
thing that catches a shared misunderstanding:

| | |
|---|---|
| Four-centre integrals and the Fock digestion | **4.3e-13** against libfint |
| Fock matrix at 2773 functions | **1e-8** against gpu4pyscf |
| Overlap, kinetic, nuclear attraction | **2.2e-15**, **1.2e-14**, **6.2e-13** against libfint |
| Two- and three-centre Coulomb | **4.5e-13** against libfint |
| Cartesian multipoles | **39 of 39** components against pyscf's `int1e_r`/`_rr`/`_rrr` |
| RHF through the public API | matches pyscf; water `-75.983974469890` |
| Folded vs enumerated digestion | **6.7e-15** on a symmetric density, where they must agree |

The multipole comparison refuses to judge a single component until the
**overlap** agrees first. That control exists because two days once went into a
disagreement that turned out to be in the oracle, not the kernel: the reference
system had a contracted shell, and pyscf renormalises contracted shells where
the reference did not. An oracle gets validated on a quantity you already trust
before it is allowed to convict anything.


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

Validation and comparison targets used during development, all open source and
all reproducible: [libcint], [pyscf], [gpu4pyscf], and [GAMESS libERI].

[libcint]: https://github.com/sunqm/libcint
[pyscf]: https://github.com/pyscf/pyscf
[gpu4pyscf]: https://github.com/pyscf/gpu4pyscf
[GAMESS libERI]: https://github.com/gms-bbg/gms_libERI

## Licence

MIT. See [LICENSE](LICENSE).
