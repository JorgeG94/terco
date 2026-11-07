# Symmetry-unique shell-pair list (sub-project)

Branch: `sym/unique-pairs`, off `symmetry` at ee2cb3d (S0-S14.4).

## The idea (from the maintainer, 2026-11-07)

> We dispatch batches of shell pairs that are assembled as quartets. We
> assemble those pairs at startup after screening. Only screened pairs are
> added to the list. Screened, symmetry-unique (point group) pairs could be
> added only, and then in the Fock the necessary missing bits can be done.

## How it differs from what is merged (S14.2, scheme 2a)

- **Today:** the pair list holds every screened pair. Bins are orbit-closed
  (`trc_bin_orbits_t`), with representatives first. A skeleton segment
  restricts one side of a bin pair to representatives (weight = orbit size)
  and keeps the other side complete. G is symmetrized once after the fold.
  The S13/S14 measurements are: water 1.47x, B2H6 3.35x, benzene 3.53x,
  naphthalene 3.79x fewer quartets. Pair data, Schwarz Q and the device
  arrays are all full size.
- **Proposed:** the stored pair list holds only the representatives, about
  1/|G| of the pairs. That shrinks the pair data, the Q and dsh tables, the
  bins and the device memory, and makes planning cheaper. The other side's
  pairs, and any other "missing bits", are produced inside the Fock build.

## What the study must settle

- **What the missing bits are.** Quartets with both sides restricted to
  representatives do not cover every quartet orbit: for a representative ab,
  the ket (g cd) is generally not a representative. The image pairs must
  therefore be reconstructed. Under an Abelian operation (axis sign flips),
  an image pair's primitive data is its representative's with the centres
  sign-flipped and the component signs from `ao_sign`. That could be done
  on the fly in the kernel's pair read (a generator change: apply the mask
  to the pair centre), or by expanding the pairs per launch on the device
  (a scratch buffer). Other routes may be better.
- The gain over scheme 2a: quartet counts (expected similar), memory,
  planning time, and the device footprint for large molecules (S14.2 noted
  that `sym_pairs_build` allocates nshell^2/2 ints).
- The interaction with screening (orbit-max Q and dsh), exchange tiles, the
  incremental Fock, the MD engines (S14.4), parity skipping
  (`sym/parity-skip`), DF and gradients.
- The compile-time cost of any generator change, and the GPU risks (register
  pressure, divergence).
- Validation: the same oracle as S14 (symmetrized G equals full G within the
  S14.2 gate), plus bit-identity with symmetry off.

## Status

Design study started 2026-11-07. It produces a design document and a CPU
prototype, then stops for maintainer sign-off.
