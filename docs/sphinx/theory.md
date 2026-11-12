# Theory

This site documents how to *use* terco: building it, driving the CLI, the
TOML input format, the output, and the Fortran/C APIs. It does not re-derive
the theory behind HF, DFT, MP2 or PCM. The exceptions are terco's own use of
point-group symmetry and the closed-shell RI-CCSD(T) implementation, each
described in its own section below.

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
projected guess, RI-MP2, RI-CCSD, CPHF/CPKS, the Hessian, Fukui indices, restarts, the
geometry optimiser and more than one MPI rank. `on`/`auto` runs these in C1;
a named group refuses them.

## RI-CCSD and RI-CCSD(T)

The closed-shell density-fitted coupled-cluster implementation
(`src/post_hf/ccsd/`) is the one method for which this page says how terco
computes it, because the choices in it (what is never stored, what is shared over
GPUs, how convergence is judged) change what a user sees. The module headers in `src/post_hf/ccsd/` are the
authority for everything below; the command-line options are in {doc}`cli`, the
input keys in {doc}`input` and the library calls in {doc}`library`.

### What it is, and what it is not

CCSD and CCSD(T) energies for a closed-shell restricted Hartree-Fock reference
with canonical orbitals, in spin-adapted (spatial-orbital) form, with the
two-electron integrals density fitted. The equations are those of PySCF's
`rccsd` and `ccsd_t` (see Validation), ported from the same author's
`mqc_czt_rcc.f90` in [metalquicha](https://github.com/JorgeG94/metalquicha).

The limits, stated plainly:

- **Closed shell only.** An unrestricted or open-shell reference is refused.
- **Hartree-Fock reference only.** A Kohn-Sham reference is refused.
- **Energies only.** There are no gradients, so no optimisation, Hessian or
  response after a converged RI-CCSD (those calls refuse), and no restart of
  the RI-CCSD itself.
- **No point-group symmetry.** A named group (`--symmetry c2v`) is refused;
  `--symmetry on` runs in C1 and says so.
- **An aufbau reference is required.** A converged SCF that holds a non-aufbau
  occupation (kept by MOM, or reached by GDM) is refused, as for RI-MP2: the
  amplitude denominators are not a ground state's.
- **No effective core potential**, and no linearly dependent functions dropped
  (the method needs every orbital of the basis).
- **Size is bounded by the memory of one GPU.** The integrals are built on
  every rank (replicated), so more GPUs speed an iteration up but do not
  make a larger system fit; see {doc}`cli` for the memory formula and its gap.
- **The speed-up from more GPUs levels off.** Only the particle-particle ladder
  and (T) are shared; W, the ring terms, the T1 and T2 equations and DIIS run
  on the first rank.

### The equations

The unknowns are `t1(i,a)` and `t2(i,j,a,b)` (`i`, `j`, `k`, `l` active
occupied, `a`, `b`, `c`, `d` virtual; chemists' notation, real orbitals, first
index fastest, which is also PySCF's order). The amplitude equations are
Eqs. (35) (T1) and (36) (T2) of S. Hirata, R. Podeszwa, M. Tobita and R. J.
Bartlett, J. Chem. Phys. **120**, 2581 (2004), with the intermediates of its
Eqs. (37)-(41), in the form PySCF's `rccsd.update_amps` (the non-CC2 branch)
and `rintermediates.py` write them. Define

```
tau(i,j,a,b) = t2(i,j,a,b) + t1(i,a) t1(j,b)
L(k,c,l,d)   = 2 (kc|ld) - (kd|lc)
```

and the one-particle intermediates `Foo(k,i)`, `Fvv(a,c)`, `Fov(k,c)`,
`Loo(k,i)`, `Lvv(a,c)` and the two-particle `Woooo(k,l,i,j)`,
`Wvoov(a,k,i,c)`, `Wvovo(a,k,c,i)`. For canonical orbitals the Fock off-diagonal
`f_ov` vanishes, the diagonal `e_p` is not carried inside the intermediates, and
the new amplitudes are divided by the orbital-energy denominators:

```
t1'(i,a)     /= e_i - e_a
t2'(i,j,a,b) /= e_i + e_j - e_a - e_b
```

The T1 equation is the sum of `Fvv t1`, `-Foo t1`, the `Fov` term, the
`[2 (kc|ia) - (ki|ac)] t1(k,c)` term, and the two `tau` contractions with the
integrals `[2 (kd|ac) - (kc|ad)]` and
`lam(k,i,l,c) = 2 (ki|lc) - (li|kc)`. The T2 equation, in
the order in which terco accumulates it, is: the base integral `(ia|jb)`;
two t1-dressed terms (A and B), each added with its transpose; the
`Woooo . tau` term; the particle-particle ladder `Wvvvv . tau`, added
last; the `Lvv` and `Loo` terms, each added with its transpose; and three ring
terms built from `2 Wvoov - Wvovo`, `Wvoov` and `Wvovo`, each one matrix
product. Every O(N^5) and O(N^6) contraction is a matrix product (GEMM) on the
device. The starting amplitudes are `t1 = 0` and `t2` = the MP2 amplitudes (so
`E_MP2` is the energy of the start), or amplitudes the caller supplies.
The correlation energy is the closed-shell functional of PySCF's `rccsd.energy`
evaluated with `L` and `tau`; terco reports its t1 t1 part and its t2 part. The T1
diagnostic is PySCF's, `sqrt(||t1||_F^2 / (2 o))`, `o` the number of correlated
occupied orbitals.

All intermediates of an iteration are built from the same input amplitudes, and
the new amplitudes go to separate arrays. A frozen core removes the first
`nfrozen` occupied orbitals from the correlated space (the default counts the
core by element); `o` below is the active count.

### Density fitting of every integral class

With the metric `J = (P|Q) = L L^T` and `B^P_pq = sum_Q (pq|Q) (L^-T)_QP`, every
integral is `(pq|rs) = sum_P B^P_pq B^P_rs`. The three-index pipeline (metric,
Cholesky factor, `(mu nu|P)` over blocks of auxiliary shells, the orbital fold
for a spherical basis) is RI-MP2's; the half transform produces the
occupied-occupied, occupied-virtual and virtual-virtual factors from one
transformed tensor. Kept on the device and the host:

- the factors `b_oo`, `b_ov`, `b_vv` (`b_vv`, of size `v^2 naux`, is the largest
  array), and the blocks `(ij|kl)`, `(ij|ka)`, `(ij|ab)`, `(ia|jb)`;
- **no** `(ab|cd)` block, ever: the ladder builds it from `b_vv` a batch of
  columns at a time;
- **no** dense `(ia|bc)` block (`o v^3` doubles): every contraction that needs it
  is done from `b_ov` and `b_vv` directly, e.g. `(kd|ac) = sum_P b_ov(kd,P) b_vv(ac,P)`
  with the `t1` contraction done first (`N1(k,a;P) = sum_d t1(k,d) b_vv(d,a;P)`).
  The reference (T) routine and the tiled one build such a copy for the duration of
  the call (`ints%fill_ovvv`).

### The particle-particle ladder

`L(i,j,a,b) = sum_cd (ac|bd) tau(i,j,c,d)` plus its t1 dressing is the largest term.
It has three forms, none holding `(ab|cd)` whole:

- **assembled**: per batch of `cd` columns, `(ac|bd)` by a strided batched GEMM
  of `b_vv` against itself, then one GEMM with `tau`;
- **direct**: per auxiliary function `P`, `L += G_P tau_ij G_P^T`, restructured into
  a batched GEMM, a permutation and one GEMM per batch of functions; it wins on
  counted multiply-adds when `v` is above `2 naux o^2/(naux + o^2)`, about `2 o^2` for
  `naux` much larger than `o^2` (not re-measured against the half form);
- **half**: `(ac|bd) = (bd|ac)`, so only the `c >= d` columns are built, in passes
  over a range of the virtual index `d`, and the `c < d` half comes from the
  symmetry of `tau`: `L(i,j,a,b) = N(i,j,a,b) + N(j,i,b,a)` when
  `tau(i,j,d,c) = tau(j,i,c,d)`, which holds for the MP2 start and every converged
  `t2`. About `naux v^4/2 + o^2 v^4/2` multiply-adds. The dressing is folded into the
  build through the dressed factor `bt(a,c;P) = b_vv(a,c;P) - sum_k t1(k,a) b_ov(k,c;P)`,
  which removes the separate `2 o^3 v^3` dressing.

`auto` takes the direct form past the count above, otherwise the half form when one
`d` of it fits its budget, else the assembled form. Over several ranks only the
half form can be cut, so `auto` is the half form and an explicit assembled or direct
form is refused.

Measured, as recorded in the module header (one GH200, water pentamer cc-pVTZ /
cc-pVTZ-RIFIT, `o` 20, `v` 265, 705 auxiliary functions, one undressed ladder per
iteration): assembled 0.336 s; half 0.150 s at `d_batch` 16. The ladder stage of the
14-iteration CCSD: 5.30 s assembled, 2.71 s half, 2.12 s half with the dressed build;
energies unchanged to the printed 1e-12.

### DIIS

DIIS (Pulay, *Chem. Phys. Lett.* **73**, 393 (1980); *J. Comput. Chem.* **3**, 556
(1982); see {doc}`citations`) acts on the concatenated `(t1, t2)`. The history holds,
for up to `diis_size` iterations (default 8; 0 or 1 switches it off), the new
amplitudes that the update produced and the step `e = new - old`. The
matrix `B_ij = <e_i, e_j>` is built one row per push through the linear-algebra
handle's dot, and the constrained system `[B -1; -1^T 0][c; l] = [0; -1]` is solved
on the host by Gaussian elimination with partial pivoting after scaling `B` to O(1).
A non-finite `B`, or a pivot under 1e-14, restarts the history (counted in the
result); a non-finite step is never pushed. DIIS starts at the first iteration.

Because `t2(i,j,a,b) = t2(j,i,b,a)`, with `p = (i,a)` and `q = (j,b)` the `t2` is a
symmetric `o v x o v` matrix. The **packed** history keeps only `p >= q`
(`o v + o v (o v + 1)/2` numbers per stored vector instead of `o v + o^2 v^2`), with
the off-diagonal elements of the error scaled by `sqrt(2)` so that the plain dot of
two packed errors is the dot of the full symmetric steps; the extrapolation writes
both `(p,q)` and `(q,p)`, so the result is exactly pair symmetric. At `o` 40,
`v` 530, size 8 (ten waters, cc-pVTZ) the module header gives 28.8 GB packed against
57.5 GB full. The **host** history keeps the vectors in host memory and maps one slot
at a time on the device. Packing is used exactly when the ladder is told the
amplitudes are pair symmetric, which is every production run.

### Convergence

The iteration stops as converged when both

```
|E_k - E_(k-1)|      < conv_energy      (default 1e-8)
max |t' - t|         < conv_amplitude   (default 1e-6, over t1 and t2)
```

hold. `E_0` is the energy of the starting amplitudes. A step or an energy that is
NaN or infinite stops the run as not converged with a message. PySCF's test is
`|dE| < 1e-7` and `normt < 1e-5`, where `normt` is the 2-norm of the difference of
the vectors `(t1, packed lower triangle of t2)`, each off-diagonal pair counted
once. terco's `step_norm2` is the 2-norm over the *full* `t2`, so with a
pair-symmetric step

```
normt <= step_norm2 <= sqrt(2) normt
```

and `step_norm2` is not PySCF's number. The tests are not nested either: `max|e| < 1e-6`
allows `normt` up to `sqrt(N') 1e-6`, `N' = o v + o v (o v + 1)/2`, which exceeds
1e-5 once `N' > 100`; terco's amplitude test is the stricter per element and can be the
looser in the 2-norm, while its energy test is the tighter. Both tolerances are
options (`--cc-econv`, `--cc-tconv`).

### The (T) correction

For each virtual triple `a >= b >= c` and each permutation `p` of `(a,b,c)`,
the (T) of Rendell, Lee and Komornicki (J. Chem. Phys. **94**, 442 (1991)), as
PySCF's `ccsd_t_slow.py` writes it, of the converged amplitudes is

```
W_p(i,j,k) = sum_f (i a_p|f b_p) t2(k,j,c_p,f) - sum_m (i a_p|j m) t2(m,k,b_p,c_p)
V_p(i,j,k) = (i a_p|j b_p) t1(k,c_p)
Z_p        = [4 X(i,j,k) + X(k,i,j) + X(j,k,i) - 2 X(k,j,i) - 2 X(i,k,j) - 2 X(j,i,k)] / D3,
             X = W_p + V_p / 2,   D3 = e_i + e_j + e_k - e_a - e_b - e_c
E(T)       = 2 sum_abc sum_p sum_ijk Z_p(i,j,k) [ sum_q W_WMAP(p,q)(perm_q(i,j,k)) ]
```

with the fitted integrals. Three routines compute it, and the energy does not depend
on the choice:

- **reference** (`--cc-triples reference`): the host implementation, the oracle the
  device routines are checked against; OpenMP over `(a,b)` pairs, deterministic. No
  memory budget.
- **tiled** (`tiled`): the device routine driven by virtual triples, with `W` tiles of
  `6 o^3` doubles per triple, sized by `--cc-triples-budget`; it runs on the host too.
- **ijk** (`ijk`): the device routine driven by occupied triples `i >= j >= k`,
  restructured so that every library call is sized by the virtual space (tall GEMMs on
  `Q`, `24 v^3` bytes per occupied triple), with a cache of integral slices. Measured
  on one GH200 (`H10O5`, cc-pVTZ / cc-pVTZ-RIFIT, `o` 20, `v` 265,
  705 auxiliary functions, default budget): the virtual-triple routine 11.90 s, the
  ijk routine 3.42 s, E(T) the same to all 12 printed digits (module header).

`auto` is the ijk routine in a GPU build and, in any build, whenever there is more than
one rank; it is the reference routine in a host build on one rank. The per-triple energies
are stored at the triple's canonical position and added in a fixed order, so E(T) does
not depend on the tile size, the slice cache or thread scheduling. Memory by routine and
the budget's meaning are in {doc}`cli`.

### Several GPUs

One MPI rank per GPU (no two ranks share a device). Every rank builds all the integrals
(replicated). Per iteration:

1. rank 0 sends the amplitudes to each other rank, point to point (the unique half of
   `t2` when pair symmetric), through device buffers when the MPI library takes them
   (`TERCO_MPI_DEVICE=0` forces the host path);
2. every rank builds `tau` and the dressed factor and runs the half-form passes of
   its range `d_lo..d_hi` of the virtual index; the ranges are cut from the counted flops
   so that every rank has the same total, rank 0 getting the smaller share (an empty one
   when its other work already exceeds an equal share, e.g. adenine cc-pVDZ at any
   rank count);
3. rank 0 meanwhile does everything that does not read the ladder: the intermediates
   and W, the T1 equation, every T2 term but the ladder, the rings;
4. one agreement on failure (any rank's refusal is every rank's, with the reason of the
   lowest-numbered rank, and no rank is left waiting);
5. each follower sends its partial ladder to rank 0, which adds it on its device, adds the
   ladder to `t2` last, applies the denominators and runs DIIS and the convergence test.

At the end rank 0's result and amplitudes are broadcast, so every rank returns the same
status and energies. The (T) is shared by occupied pairs `(i,j)`, `i >= j`, dealt by
decreasing cost to the part with the fewest triples; each rank computes its triples,
the per-triple vectors are summed with `MPI_SUM` (each entry has one non-zero term, so
the sum is exact) and E(T) is summed in canonical order on every rank. Only the ijk
routine does this; an explicit `reference` or `tiled` is refused over several ranks.

**Order of summation.** The ladder is summed in a different order than one rank sums
it, and is added to `t2` after the rings in every run (one rank included), so the
RI-CCSD energy over several ranks agrees with the one-rank value to rounding, not
bitwise, and the iteration count is the same. The multi-rank SCF in front of it
differs from a one-rank SCF at about 1e-14 Eh. `test/rccsd_mpi.F90` bounds the
differences at 1e-14 for `E_CCSD`, every amplitude and `E(T)`; measured on water and
hydrogen peroxide, host and GPU builds, 2 to 4 ranks: 1.7e-16 at most for `E_CCSD`,
4.2e-17 for the amplitudes, 3.5e-18 for `E(T)`.

Measured iteration times are in {doc}`cli` (Several GPUs): the CCSD iteration of adenine
cc-pVDZ takes 0.183 s on one V100 and 0.119 s on two; the water hexamer cc-pVTZ is 2.87
times faster on four GPUs than on one. The gain levels off with the rank count because
W, the ring terms, the T1 and T2 equations and DIIS stay on rank 0.

### Validation

The energies are checked against PySCF's DF-RCCSD(T) (`pyscf` never runs in the test
suite; its numbers are frozen in `test/pyscf_refs/rccsd_ints_*.ref`, written by
`generate.py`), on water and a C1 hydrogen peroxide, cc-pVDZ with cc-pVDZ-RIFIT, spherical,
frozen core 1 and 2. With the references converged to `conv_tol` 1e-12 and normt 1e-10 and
terco converged to 1e-12 / 1e-10 (SCF 1e-13, commutator 1e-11), the tests require:

| Quantity | Tolerance (Eh) | Test |
|---|---|---|
| E_MP2 | 1e-10 | `check_capi_rccsd` |
| E_CCSD correlation | 1e-10 | `check_rccsd_run`, `check_capi_rccsd`, `check_cli_rccsd` |
| E(T) | 3e-12 | `check_rccsd_triples`, `check_capi_rccsd`, `check_cli_rccsd` |
| CCSD(T) total energy | 3e-11 | the same |

Other checks, each against something that can fail: the fitted integral blocks as
sign-invariant probes (norms and contractions with the MO dipole integrals); one
amplitude update from fixed amplitudes against PySCF's `update_amps` and against an
independent host implementation written term by term from the PySCF einsums; the (T) at
fixed amplitudes (independent of the CCSD solver) and against an independent evaluation
over all `(a,b,c)` with no `a >= b >= c` restriction; the device (T) routines against the
host reference to 1e-16 and for tile independence (forced tiles of 1, 7 and 37 triples);
DIIS on against off, assembled against direct ladder; and the multi-rank run against the
one-rank run (above). The command-line document states the agreement a user sees: within
5e-12 Eh for the CCSD(T) total at tight thresholds, about 1e-8 Eh at the default ones.

### Not in this implementation

Nothing here is a gradient, an open-shell or Kohn-Sham reference, a symmetry-adapted
calculation, a CC2 branch, or a level shift of the denominators (the
intermediates leave the Fock diagonal out and the update subtracts the same diagonal, so
a shift would have to enter both).
