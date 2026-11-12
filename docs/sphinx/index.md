# terco

A GPU electron-repulsion-integral and self-consistent-field library,
written in standard Fortran. Compute kernels are `do concurrent`; OpenACC
moves data and nothing else. terco's claim is that standard Fortran is
competitive with hand-written CUDA for this kind of code.

This is the **user documentation** for terco v0.1.0: how to build it, how to
drive the `terco` command-line program, the TOML input file format, what the
output looks like, and how to call terco as a library (Fortran API or C
ABI). It does not re-derive the theory behind HF/DFT/MP2/PCM; see
{doc}`theory` for that. {doc}`theory` also describes how terco implements RI-CCSD(T).

```{toctree}
:maxdepth: 2
:caption: Contents

overview
install
cli
input
output
library
api
functionals
theory
citations
ai_policy
```
