# Theory

This site documents how to *use* terco: building it, driving the CLI, the
TOML input format, the output, and the Fortran/C APIs. It does not re-derive
the theory behind HF, DFT, MP2 or PCM. The one method it does describe is
terco's own use of point-group symmetry, below.

For the underlying theory, the methods terco shares with its sister project
metalquicha are documented in metalquicha's own user guide:

<https://metalquicha.readthedocs.io/en/latest/>

terco's own method-level references -- the papers the ERI algorithms and
screening/binning design are drawn from -- are listed in the top-level
[README](https://github.com/JorgeG94/terco#references), with BibTeX
entries for all of them.

## Point-group symmetry

Symmetry is opt-in (`--symmetry`, see {doc}`cli`). terco uses it three ways:
the SCF is solved block by block in symmetry-adapted functions, the Fock
matrix and the four-centre gradient are built from a reduced set of integral
quartets, and integral components that vanish by symmetry are not computed.
None of these changes the answer beyond round-off.

### Groups

terco works with the Abelian groups D2h, D2, C2v, C2h, Cs, C2, Ci and C1.
Every operation is a diagonal sign matrix in the group's standard frame, so
it permutes the atoms and maps each Cartesian AO onto another AO times a
sign. A molecule of a higher group (a symmetric or spherical top, a linear
molecule, an atom) is run in its largest Abelian subgroup and frame, the ones
pyscf uses. A group named by the user must be the detected group or one of
its subgroups; D2h to D2 is not offered, as pyscf does not offer it.

### Detection and symmetrisation

The geometry is shifted to the centre of nuclear charge and the
charge-weighted second-moment tensor is diagonalised. When its eigenvalues
are distinct, the seven sign masks are tested in the principal frame: a mask
is a symmetry when each atom maps onto an equivalent atom (same element, same
ghost status and same basis) within `--symmetry-tol` (1e-5 Bohr). The
passing masks form the detected group. The molecule is then symmetrised,
atoms on a symmetry element getting exact zero coordinates, and the largest
atom displacement this made is reported. With `on`/`auto` the results are
rotated back so that they come out in the input frame; with a named group they
stay in the standard frame.

### The SO-blocked SCF

The AO basis is split into symmetry-adapted functions (SOs), one block per
irrep. Each block is orthogonalised on its own (canonical orthogonalisation,
with the same normalisation and linear-dependence threshold as the C1 code)
and diagonalised on its own, so the eigenproblem is several small ones
instead of one large one. The orbitals of all blocks are merged by energy, so
the aufbau occupation is the first `nocc` columns, as in C1, and every
orbital carries its irrep label. `--irrep-nelec` instead fills the lowest
orbitals of each named irrep. The density and the Fock matrix are
symmetrised each iteration.

### The skeleton Fock build

The Fock build is the expensive part, and it is where symmetry saves most.
Let the group act on a shell quartet $q$. For a density $D$ invariant under
the group, the contribution of the image quartet is the contribution of $q$
transformed by the AO signed permutation $R(g)$:
$\mathrm{Dig}_{gq}(D) = R(g)\,\mathrm{Dig}_q(D)\,R(g)^T$. Take one
representative $q$ of each orbit of quartets, weight it by the size of its
orbit $w_q = |G|/|\mathrm{Stab}(q)|$, and sum:
$G_\mathrm{skel} = \sum_{q\in Q_0} w_q\,\mathrm{Dig}_q$. Then

$$
G(D) = \frac{1}{|G|}\sum_{g} R(g)\,G_\mathrm{skel}\,R(g)^T ,
$$

because the sum over the images of $q$ is $|\mathrm{Stab}(q)|^{-1}$ of the sum
over the whole group, which the weight undoes. The skeleton matrix itself is
not the Fock matrix (it differs by order one), and the symmetrisation, the same
one that is applied to the density, is what makes it exact. The identity
holds for J and K separately, for the density-weighted Schwarz screen (the
screening bounds are taken as the maximum over each pair orbit, so a quartet
and its images pass or fail together), for exchange-tile screening and for
the incremental Fock build (the density difference is invariant too).
Measured, the symmetrised skeleton matrix agrees with the full one to 3e-14
(water) up to 2e-12 (naphthalene, cc-pVDZ) on elements as large as 24: the
summation noise of 2e7 quartets.

In the engine the bra pairs of each bin pair are cut to orbit
representatives and the ket runs over the whole bin, with the orbit weight
carried on the launch.

What it saves is the orbit reduction of the quartet count, at most the
effective order of the group. An operation that fixes every atom (the
molecular plane of a planar molecule) acts on the shells as a pure sign and
reduces nothing here, so planar D2h molecules reach 4 of 8 and only
nonplanar ones such as B2H6 reach 8. Measured quartet reduction (cc-pVDZ,
pair screening 1e-10, density-screened):

| System | Group | atoms / atom orbits | Quartets, full over skeleton |
|---|---|---|---|
| water | C2v | 1.50 | 1.47x |
| SF6 | D2h | 1.75 | 2.44x |
| B2H6 | D2h | 2.67 | 3.34x |
| benzene | D2h | 3.00 | 3.53x |
| naphthalene | D2h | 3.60 | 3.79x |
| pentacene | D2h | 3.60 | 3.85x |

Density fitting and several MPI ranks take the full build.

### Parity skipping

Where all four centres of a quartet lie in a coordinate plane (or on an
axis) of the molecule, every component integral of odd combined parity in
that coordinate is exactly zero. The group tells which axis can be used:
one whose mirror the group contains and in which atoms lie exactly. The
engine permutes the molecule's axes so that this one is the kernels' z axis,
bins the shells in the plane apart, and runs masked kernel variants for the
classes that carry the cost, which leave the zero components out from the
primitive recurrences on. By operation count this removes about 30 percent of
the work of a planar D2h molecule (34 to 52 percent of the component
integrals are zero), on top of the orbit reduction; it is a model, not a
timing. The zeros are exact, so the result is unchanged. `TRC_PARITY=0` in the
environment turns it off.

### The skeleton gradient

The two-electron gradient is a sum over quartets of a force on the (at most
four) centres of each. Rotating a quartet by $g$ rotates its force and moves
it to the image atoms, so the same representatives and weights give
$\nabla_\mathrm{skel}=\sum_{q\in Q_0}w_q F_q$, and the gradient is its
projection onto the symmetric forces,
$\nabla = |G|^{-1}\sum_g P_g\,\nabla_\mathrm{skel}$, where $P_g$ permutes the
atoms and flips the signs of the components. The unprojected skeleton
gradient is wrong by about the size of the gradient itself; the projection
is applied to the four-centre term before it is added to the rest. Measured
against the full gradient it agrees to 3e-15 (water) to 2e-12 (benzene) on
random densities, and end to end to 4e-15..1e-14. On a V100 the gradient's
integral kernels fell by 1.56x (benzene), 2.60x (naphthalene) and 3.31x
(pentacene), less than the quartet reduction because the skeleton launches
are more numerous and shorter.

### The size rule

A skeleton build computes one quartet per orbit and weights it, so it saves
the orbit reduction and pays for extra, smaller launches. The saving depends
on how large the atom orbits are, not on the basis or the atom count alone: a
molecule whose atoms sit on many symmetry elements at once has small orbits
(SF6's sulfur is fixed by all eight operations). With neither
`--sym-skeleton` nor `--no-sym-skeleton`, terco uses the skeleton build,
for the Fock matrix and the gradient together, when the number of atoms is at
least twice the number of atom orbits. Measured on a V100 (cc-pVDZ,
`--thresh 1e-12`, ERI time per Fock build against `--symmetry off`):

| System | atoms | atom orbits | atoms / orbits | ERI time, symmetric over off | Build chosen |
|---|---|---|---|---|---|
| SF6 | 7 | 4 | 1.75 | 0.44x (slower) | full |
| benzene | 12 | 4 | 3.00 | 1.41x | skeleton |
| naphthalene | 18 | 5 | 3.60 | 2.07x | skeleton |
| pentacene | 36 | 10 | 3.60 | 2.78x | skeleton |

The threshold is set at the low end of the interval the data leaves open
(between 1.75 and 3.00), so only what was measured to lose is declined. It
is provisional: B2H6, at 2.67, lies in the interval, and the rule has not been
re-timed on a device. Groups of order 2 and 4 are declined wherever an atom
sits on a symmetry element, which is most molecules (water 1.47x); the
`--sym-skeleton` and `--no-sym-skeleton` flags, or `TRC_SYM_SKELETON=0/1`,
override it.

### Accuracy

The symmetric and the C1 SCF solve the same problem on the same symmetrised
geometry, and the total energies agree to 3e-13..1.4e-12 Eh against a 1e-10
convergence gate, with identical iteration counts. The gradient agrees to the
same round-off. The exception is DFT with `on`/`auto`: the grid is built on
your axes and the Fock matrix symmetrised afterwards, so energies differ from
`--symmetry off` at second order in the grid error and gradients at first
order (1e-7 to 1e-6 on water at the default grid level 3, 2e-9 at level 5; see
{doc}`input`). A named group puts the grid on the group's axes, so a DFT energy
on a rotated input moves with the grid orientation.

### Not supported with symmetry

Density fitting, PCM, fractional occupations, SOSCF, MOM, GDM, a given or
projected guess, RI-MP2, CPHF/CPKS, the Hessian, Fukui indices, restarts, the
geometry optimiser and more than one MPI rank. `on`/`auto` runs these in C1;
a named group refuses them.
