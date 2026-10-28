# MD engine benchmarks (w5 / w20 cc-pVQZ)

Harness for terco's McMurchie-Davidson J/K engine against GPU4PySCF, as
measured on a V100 (2026-10-27/28) and to be repeated on GH200.

| file | what |
|---|---|
| `w5.xyz`, `w20.xyz` | water clusters (15 / 60 atoms) |
| `cc-pvqz-pyscf17.json` | pyscf 2.14's cc-pVQZ for H through Cl at 17 digits, so both codes see the same basis |
| `run_terco.sh` | terco iteration-4 energy and Fock time; `nsys` to profile |
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

## On another GPU

The cuBLAS GEMM path for MD exchange was measured on three GPUs and lost on
all of them (w20 cc-pVQZ per Fock build: V100 26.4 vs 20.5 s, A100 8.23 vs
7.76 s, GH200 4.30 vs 4.07 s, FP64 tensor cores notwithstanding), so it was
removed; MD exchange is the fused `do concurrent` path only.

1. Build: `-DTERCO_GPU_ARCH=<ccXX> -DTERCO_LMAX=4 -DTERCO_NVTX=ON` (compile
   for the GPU you run on: a cc70 build on a Hopper runs, slowly).
2. `B=<build> GPU=0 ./run_terco.sh`, then `./run_terco.sh nsys`.
3. `python g4p_jk.py w20.xyz cc-pvqz 1e-12 > g4p_w20.log`, then
   `python join_profile.py w20.sqlite g4p_w20.log` (or with no log for
   terco's per-class times alone).
4. The Rys/MD crossover (`TERCO_RYS_G_LTOT`, 9 on V100) and the chunk plans
   (`scripts/gen_rys.py` FOCK_CHUNK_PLAN, `--chunk-*` flags) were raced on
   a V100 and held on cc-pVQZ; re-race them before tuning for a new GPU.
5. Multi-GPU: `mpirun -np N` on the same inputs (every engine splits over
   ranks; energies must not change with N).
