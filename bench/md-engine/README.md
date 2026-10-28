# MD engine benchmarks (w5 / w20 cc-pVQZ)

Harness for terco's McMurchie-Davidson J/K engine against GPU4PySCF, as
measured on a V100 (2026-10-27/28) and to be repeated on GH200.

| file | what |
|---|---|
| `w5.xyz`, `w20.xyz` | water clusters (15 / 60 atoms) |
| `cc-pvqz-pyscf17.json` | pyscf 2.14's cc-pVQZ for H through Cl at 17 digits, so both codes see the same basis |
| `run_terco.sh` | terco iteration-4 energy and Fock time, MD exchange `fused` vs `gemm`; `nsys` to profile |
| `g4p_jk.py` | GPU4PySCF `get_jk` at the converged density plus its per-class timer |
| `join_profile.py` | per-class terco vs GPU4PySCF table from the two |

Run at the same screening threshold: terco's default `--thresh 1e-12`, and
`python g4p_jk.py w20.xyz cc-pvqz 1e-12` (GPU4PySCF's default is 1e-13).
pyscf uses CODATA 2010 Bohr and terco 2018; the w20 energies agree to ~1e-10.

Basis DATA matters as much as the basis name. pyscf's cc-pVQZ has the
redundant primitives of its general contractions removed; the BSE files
(terco's and metalquicha's `cc-pvqz.json`) keep them. Same 2300 functions
and 700 shells for w20, but 1140 primitives (at most 9 per shell) against
1440 (at most 12), about 1.5x the Fock time, for the same energy (4e-11).
Compare codes on this file; time terco's defaults on the BSE ones.

## V100 reference (s per full J+K build, RHF, 1e-12)

| system | terco (f364ff8) | GPU4PySCF 1.8.1 |
|---|---|---|
| w5 cc-pVQZ | 0.77 | 1.11 |
| w20 cc-pVQZ | 13.16 | 11.6 |
| w20 cc-pVQZ, A100 (GPU4PySCF 1.9.0) | 7.76 | 9.22 |

Energies (terco iteration 4, `--iterations 4`): w5 -380.228026426339,
w20 -1520.929072580745. Converged w5 RHF: -380.228710965728.

## GH200 plan

The V100 has no FP64 tensor cores, so the fused (`do concurrent`) MD
exchange beat the cuBLAS GEMM path (w20: 20.5 vs 26.4 s) and Rys beat MD
for g classes up to total L 9. GH200 runs DGEMM on FP64 tensor cores
(~67 vs ~34 TF/s vector), so both choices need re-measuring.

1. Build: `-DTERCO_GPU_ARCH=cc90 -DTERCO_LMAX=4 -DTERCO_NVTX=ON`, plus
   `-DTERCO_MDG_CLASSES=<d/f tags>` to also generate GEMM kernels for d/f
   classes (route them at run time with `TRC_MDG_CLASSES=...`).
2. `B=<build> GPU=0 ./run_terco.sh` (fused vs gemm), then with `nsys`.
3. `python g4p_jk.py w20.xyz cc-pvqz 1e-12 > g4p_w20.log`, then
   `python join_profile.py fused_w20.sqlite g4p_w20.log`.
4. If GEMM wins: re-race `TERCO_RYS_G_LTOT` (9 on V100) and the chunk
   plans (`scripts/gen_rys.py` FOCK_CHUNK_PLAN, `--chunk-*` flags).
5. Multi-GPU: `mpirun -np N` on the same inputs (every engine splits
   over ranks; energies must not change with N).
