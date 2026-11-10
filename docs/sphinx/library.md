# Using terco as a library

terco is a library first; the `terco` command-line program (see {doc}`cli`)
is built on top of the same Fortran API described here. See {doc}`api` for
the generated reference of every public entry point.

## From Fortran

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

The basis, the pair list and the ERI object are built **once per geometry**
and reused by every Fock build; rebuilding them per SCF iteration throws
away the container's entire reason for existing.

`fock_many` takes `(ndens, nao, nao)` -- the density index **first**: with
it last, the N atomic updates for one `(mu,nu)` sit `nao^2` apart and the
batch stops paying past N = 4.

## From another compiler (the C ABI)

A gfortran or ifx host cannot read nvfortran's `.mod` files, so it goes
through the C ABI. `libterco.so` is the default build -- nvfortran links
it, so it already carries the NVHPC and OpenACC runtimes as its own
`DT_NEEDED`, and the host's link line needs nothing but `-lterco`.

From C, `#include "terco.h"` (generated from `trc_capi.F90`, under
`include/`) and link against `libterco.so`:

```c
#include "terco.h"

trc_handle_t h;
double z[3]   = {8.0, 1.0, 1.0};          /* O, H, H */
double xyz[9] = { ... };                  /* Bohr, atom-major x,y,z */
double energy;

trc_create(&h);
trc_set_molecule(h, 3, z, xyz, /*charge=*/0, /*multiplicity=*/1);
trc_set_basis_json(h, "6-31g.json");
trc_set_method(h, "", /*grid_level=*/3);  /* empty name: Hartree-Fock */
trc_set_convergence(h, 1.0e-10, 1.0e-7, /*max_iter=*/100, /*diis=*/10);
trc_run_scf(h);
trc_energy(h, &energy);
trc_destroy(h);
```

This is the same sequence `test/link_c.c` runs end to end (water/6-31G RHF)
as part of terco's own test suite, against `terco.h` and `libterco.so`
built by the project's own C compiler -- not nvfortran, not gfortran -- so a
mismatch between the generated header and the real ABI shows up there
before any external caller sees it.

From a Fortran caller built with a *different* compiler than the one that
built `libterco.so`, compile `include/trc_c_interfaces.f90` with your own
compiler and `use trc_c_interfaces`; that file is the authoritative
declaration of the ABI and should not be re-typed by hand:

```fortran
use trc_c_interfaces
rc = trc_basis_create_libcint(atm, natm, bas, nbas, env, nenv, 1, hbas)
rc = trc_eri_create(hbas, 1.0e-10_dp, heri)
rc = trc_fock(heri, hbas, d, g, 1.0d0, 1.0d0, 1, 0.0d0)   ! last: omega, 0 = Coulomb
```
