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

## Point-group symmetry (C ABI)

Symmetry is off until a caller asks (ABI 1.23). Set the request before
`trc_run_scf`; read the results after it.

```c
trc_set_symmetry(h, "on", 0.0);        /* tol 0: the default, 1e-5 Bohr */
trc_run_scf(h);

char group[16];
trc_symmetry_group(h, group, sizeof group);     /* "C2v" */
int nirrep;  trc_symmetry_nirrep(h, &nirrep);
int *mo_irrep = malloc(nao * nspin * sizeof(int));
trc_mo_irreps(h, mo_irrep);                     /* 1-based irrep per orbital */
```

| Function | Version | Meaning |
|---|---|---|
| `trc_set_symmetry(h, group, tol)` | 1.23 | the request: `"off"` (the default), `"on"`/`"auto"` (the largest Abelian group; unsupported paths run in C1), or a group name (strict). `tol` in Bohr, 0 for the default |
| `trc_set_irrep_nelec(h, spec)` | 1.23 | electrons per irrep, `"A'=6,A\"=2"` (restricted) or `"A'=4/3,A\"=1/1"` (alpha/beta); `""` restores aufbau. Makes the request strict |
| `trc_set_symmetry_report_frame(h, frame)` | 1.23 | `1` input frame, `0` standard frame, `-1` the default (below) |
| `trc_set_symmetry_skeleton(h, on)` | 1.25 | `1` force the skeleton Fock and gradient builds, `0` force the full builds, `-1` let the size rule choose (the default). Any other value is `TRC_ERR_BADARG`. The environment variable `TRC_SYM_SKELETON=0/1` overrides all three |
| `trc_symmetry_group(h, buf, n)` | 1.23 | the group the last run used (`"C1"` when none) |
| `trc_symmetry_detected(h, buf, n)` | 1.23 | the largest Abelian group found in the geometry |
| `trc_symmetry_residual(h, &res)` | 1.23 | the largest atom displacement symmetrisation made, Bohr |
| `trc_symmetry_frame(h, rot, shift)` | 1.23 | `x_std = rot (x - shift)`: `rot` is 9 doubles row-major, `shift` the centre of nuclear charge (identity and zero for C1) |
| `trc_symmetry_report_frame(h, &frame)` | 1.23 | the frame the getters return now (`1` input, `0` standard) |
| `trc_symmetry_nirrep(h, &n)` | 1.23 | number of irreps of the group used |
| `trc_symmetry_irrep_name(h, i, buf, n)` | 1.23 | name of irrep `i` (1-based) |
| `trc_symmetry_irrep_dim(h, dim)` | 1.23 | symmetry-adapted functions per irrep (`nirrep` ints, summing to `trc_nao`) |
| `trc_mo_irreps(h, irreps)` | 1.23 | irrep of each orbital, 1-based, laid out like `trc_mo_energies` (`nao` by `nspin`, column-major); 0 for an orbital the linear-dependence cut dropped |
| `trc_geometry(h, xyz)` | 1.23 | the geometry the context holds now (Bohr, atom-major): after a symmetric run, the symmetrised one in the reported frame |

The getters return `TRC_ERR_STATE` unless the last run was symmetric (the
group-name getters unless a symmetry request ran). A request longer than 32
characters, or an irrep specification longer than 256, is refused.

**Frame.** With `"on"`/`"auto"` the geometry, gradient, orbitals, density,
multipoles and the borrowed basis and ERI handles are reported in your
**input frame** (the input symmetrised in it), rotated back exactly from the
frame the SCF ran in. With a named group they are in the **standard frame**
of the group (centre of nuclear charge at the origin, axes along the
group's): `trc_symmetry_frame` maps your input to it.
`trc_set_symmetry_report_frame` overrides either default. Energies, orbital
energies, irrep labels, occupations and atomic charges do not depend on the
frame.

**Unsupported paths.** An `"on"`/`"auto"` request runs these in C1 and logs
one line naming the reason; a named group (or any request with
`trc_set_irrep_nelec`) refuses them at the run with `TRC_ERR_UNSUPPORTED`:
density fitting, PCM, fractional occupations, SOSCF, MOM, GDM, a given or
projected guess, more than one rank, RI-MP2, CPHF/CPKS, the Hessian, Fukui
indices, restarts and the geometry optimiser.
