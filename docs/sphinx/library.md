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

## RI-CCSD and RI-CCSD(T) through the C ABI

Closed-shell RI-CCSD, and (T) on request, are in the C ABI since version 1.28
(`TRC_ABI_MINOR`); the Fortran facade `trc_api` does not carry them. They run on
the converged restricted Hartree-Fock reference in the handle, with the auxiliary
basis of `trc_set_aux_json` (the one RI-MP2 reads). Six entry points, documented in
`include/terco.h` (and in {doc}`api`):

| Call | What it does |
|---|---|
| `trc_set_rccsd(h, triples, nfrozen, aux_block)` | `triples` 1 asks for (T); `nfrozen` frozen core orbitals, -1 counts the core by element; `aux_block` auxiliary functions per three-index block (0: all at once; the energy does not depend on it) |
| `trc_set_rccsd_convergence(h, conv_energy, conv_amplitude, max_iterations, diis_size)` | the stopping rule; each 0 keeps the default (1e-8, 1e-6, 100, 8); both thresholds must hold |
| `trc_set_rccsd_algorithm(h, ladder, triples_routine, triples_budget_mib)` | the ladder form, the (T) routine and the memory of its working tile; none changes the energy |
| `trc_run_rccsd(h)` | the CCSD and, if asked, the (T) |
| `trc_rccsd_energy(h, ...)` | `e_mp2`, `e_ccsd`, `e_t`, the T1 diagnostic, the iteration count, `converged` and `has_triples`, all energies being correlation energies |
| `trc_rccsd_dimensions(h, ...)` | the frozen, active occupied, virtual and auxiliary counts of the last run |

```c
double e_mp2, e_ccsd, e_t, t1_diag, e_total;
int iterations, converged, has_triples;

trc_set_aux_json(h, "cc-pvdz-rifit.json");
trc_set_convergence(h, 1.0e-10, 1.0e-8, 100, 10);   /* the SCF */
trc_run_scf(h);

trc_set_rccsd(h, /*triples=*/1, /*nfrozen=*/-1, /*aux_block=*/0);
trc_set_rccsd_convergence(h, 1.0e-10, 1.0e-8, 0, 0); /* 0 keeps iterations, DIIS */
int rc = trc_run_rccsd(h);                           /* TRC_OK, or see below */

trc_rccsd_energy(h, &e_mp2, &e_ccsd, &e_t, &t1_diag,
                 &iterations, &converged, &has_triples);
trc_energy(h, &e_total);    /* after a converged run: the CCSD(T) total */
```

After a converged run `trc_energy` returns the CCSD (or, with (T), the CCSD(T))
total energy, not the SCF energy; RI-MP2 and RI-CCSD share it, and the later run
wins (a successful `trc_run_rccsd` drops a stored RI-MP2 result and the reverse). A
refused run leaves the other method's result alone and drops its own method's
earlier one.

What a call returns:

- `TRC_ERR_STATE`: no SCF has run, it did not converge, or no auxiliary basis was set.
- `TRC_ERR_UNSUPPORTED`, with the reason in `trc_message`: an unrestricted or
  open-shell reference, a Kohn-Sham one, an effective core potential, dropped
  linearly dependent functions, a non-aufbau reference (mom or gdm), a named point
  group, a basis above this build's density-fitting ceilings, or, over several ranks,
  an explicit reference or tiled (T) routine. `TRC_ERR_BADARG` when the frozen core
  leaves no occupied orbital.
- `TRC_ERR_NOCONV` when the iteration does not converge; the last energies are in
  `trc_rccsd_energy` with `converged` 0, (T) is not run and `trc_energy` is the SCF
  energy.
- When the CCSD converged and (T) then failed (for example for device memory, or a
  budget below one triple) the return is `TRC_ERR_UNSUPPORTED` with the (T) message;
  `trc_rccsd_energy` still gives the converged CCSD with `has_triples` 0 and
  `trc_energy` is the CCSD total.

After a converged run `trc_run_gradient`, `trc_run_hessian` and
`trc_run_polarizability` refuse. Over MPI ranks (`trc_set_comm`, one rank per device,
every rank calling `trc_run_rccsd` together) the ladder of every iteration and the (T)
are shared and every rank returns the same status and energies; see {doc}`theory`.

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
projected guess, more than one rank, RI-MP2, RI-CCSD, CPHF/CPKS, the Hessian, Fukui
indices, restarts and the geometry optimiser.
