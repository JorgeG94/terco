# Skipping integrals that vanish by symmetry inside a shell (parity skipping)

Branch: `sym/parity-skip`, off `symmetry` at b49215b. It is a separate project
from S14 by maintainer decision 1 in `docs/SYMMETRY_FOCK.md` §6a.

## Why

The S14 skeleton Fock build (scheme 2a) skips work at the **shell** level:
it evaluates one representative of each orbit of shell pairs under the
point group. An operation that maps every shell to itself gives no shell-level
saving. The clearest case is the molecular-plane mirror of a planar molecule.
At the shell level such an operation acts only as a sign on the Cartesian (or
spherical) components, (-1)^(component parity). So the gain on planar
molecules stops at |G|/2: about 3.5-3.9x of an ideal 8x for planar D2h at
cc-pVDZ (`docs/SYMMETRY_FOCK.md` §2.2).

For such an operation, an integral (ab|cd) whose combined component parity
is odd vanishes **exactly**. Today the generated kernels compute those
integrals and then discard them (or add zeros). Skipping them would recover
the missing factor, up to about 2x on planar D2h.

## Scope (to be designed; nothing decided beyond this)

- Where the parity information enters: per quartet class, per segment, or
  per launch, and whether it changes the generated per-class kernels
  (`scripts/gen_*.py`) or only their dispatch.
- Which engines are affected: HGP/Rys four-centre, the MD J and exchange
  engines, DF three-centre.
- Compile-time cost of any generator change (the per-class kernels already
  take 15-20 minutes on GPU) and register pressure.
- Interaction with S14 (scheme 2a): parity skipping should compose with
  orbit-representative bins, with the weights unchanged.
- Validation: the same oracle as S14 (symmetrized G equals full G), plus
  bit-identity with symmetry off.

## Status

Not started. The first step is a design spike in the style of S13: measure
how many evaluated integrals are parity-zero on planar test molecules,
prototype on the CPU, and recommend a scheme for maintainer sign-off.
