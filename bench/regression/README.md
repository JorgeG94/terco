# Performance-regression bench

`run.py` runs a handful of fixed inputs through a built `terco`, takes the
median stage times over a few repeats and compares them with a golden CSV
recorded earlier **on the same architecture**. It answers one question: did
this change make terco slower, or change an energy, on this GPU?

It is not:

- the cross-code comparison (terco against GPU4PySCF and TeraChem): that lives
  in the separate `terco-fock-bench` harness and uses different settings;
- a CI job. Runners have no GPU. You run it by hand (or an agent does) before
  a performance-sensitive PR and paste the table into the PR body.

Python 3.11 or newer, standard library only.

## Running it

GPU build (a device `terco`; choose the GPU explicitly, **never GPU 0 by
default**: `ctest` pins GPU 0 and other work runs there):

```sh
python3 bench/regression/run.py --terco build/app/terco --gpu 2
python3 bench/regression/run.py --terco build/app/terco --gpu 2 --quick --repeats 1   # fast check
```

`--gpu N` sets `CUDA_VISIBLE_DEVICES=N` for the runs; an already-set
`CUDA_VISIBLE_DEVICES` works instead. Without either the script refuses to
start. `--gpu 0` needs `--allow-gpu0`.

CPU node (a host build; no GPU option needed, its own golden `cpu-<model>.csv`):

```sh
python3 bench/regression/run.py --terco build-cpu/app/terco --cases fock_dz --repeats 1
```

Other options: `--cases a,b`, `--list`, `--dry-run` (print the commands),
`--repeats N` (default 3), `--warmup N`, `--out DIR` (keep the raw JSON and log
of every run plus `result.csv` in the golden's format), `--goldens DIR`,
`--arch SLUG`, `--allow-build-mismatch`, `--self-test`. Relative `--out`,
`--goldens`, `--config` and `--terco` paths are made absolute first (the runs
happen inside the output directory). Without `--terco` the newest
`build*/app/terco` in the repo is used. Cases the build cannot run (`needs` in
`cases.toml`: `lmax3` for f functions) are skipped with a note.

### What is judged and what is only reported

A case lists, under `judge` in `cases.toml`, the stages that can fail: the
per-stage means and the device-dominated stage totals (`scf.fock_mean_after_first`,
`dft.xc`, `dft.vv10`, the DF-JK and RI-MP2 steps, `pcm.per_iteration`,
`gradient.total`, ...), chosen so that one regression gives one FAIL row. Every
other recorded stage (`total_wall`, `scf.wall`, `scf.setup`, `scf.fock_first`,
`dfjk.total`, which equals the Fock total, ...) is reported with verdict
`INFO` and never fails: they hold host work, one-off device setup and file
I/O, and flake. `--list` shows the judged stages of each case.

The table has one row per recorded stage, plus rows for the energy, the
gradient RMS of gradient cases and the SCF iteration count:

- the **energy** (and the gradient RMS): **every repeat** is compared with the
  golden, and a repeat beyond `energy_tol` (1e-8 Eh; `gradient_tol` 1e-6) is a
  FAIL whatever the timings say, so a faster wrong answer fails. Repeats that
  disagree with each other beyond the tolerance FAIL as well;
- the **iteration count** against the golden is a WARN of its own (a FAIL for
  a case with `iteration_verdict = "FAIL"`, the default-incremental case, where
  the schedule of full and incremental builds is what is guarded). Repeats that
  took different iteration counts get the same verdict: an SCF that converges
  near its threshold takes n or n + 1 iterations from run to run on a GPU
  (measured: gly100 in 6-31G, 15 or 16), and the per-build means do not depend
  on the count.

Verdicts: `OK`; `WARN` (slower than the warn ratio, or a stage or case in the
golden that is missing now); `FAIL` (slower than the fail ratio);
`FAST` (more than the fail ratio faster: information, with a hint to
re-record the golden in a PR that states the intended speed-up, because a
later regression hides under a stale golden); `SMALL` (golden **and** new time
both below `min_seconds`, 0.05 s: a stage that grows from 0.04 s to 2 s is
judged, not hidden); `INFO`; `NEW` (a stage or case in the run but not in the
golden); `no golden`; and for a case that did not complete `ERROR`,
`NOT CONVERGED` or `BUILD MISMATCH`.

**Thresholds.** The configured ratios are `fail` 1.10 and `warn` 1.05. A golden
row recorded with at least 5 repeats also carries its measured spread, and its
fail ratio is then `max(configured fail, 1 + 4 x (max - min) / median)`, the warn
ratio half-way between 1 and that. With fewer repeats the configured ratios
apply. The `warn/fail` column prints the pair used for each judged stage, with a
`*` when it came from the spread.

**Errors.** One case that fails, times out (the whole process group is
killed, so no GPU process is left behind), does not converge (terco exits 2
with its JSON written) or writes no JSON does not stop the bench: it gets an
`ERROR` / `NOT CONVERGED` row naming its saved log, the other cases still run
and are judged, the table is printed, and the exit status is 2.

**Build mismatch.** A golden is only meaningful against the build it was
recorded with. Per case, if the golden's device, GPU architecture flag, `LMAX`,
`DF_LMAX`, `DF_LMAX_AUX`, Cartesian order, MPI or multicore setting, or the
compiler family differ from the binary's `terco --version`, the case is
refused (`BUILD MISMATCH`, saying exactly what differs, exit 2), and `--record`
refuses to write into it, unless you pass `--allow-build-mismatch`.

**Warm-up.** `--warmup N` runs each case N extra times first and discards them
(default 1 for a GPU build, 0 for a host build): the first run on a device pays
module loading and kernel setup.

Exit status: 0 no regression, 1 at least one FAIL, 2 usage error or at least
one case that errored, did not converge or was refused (2 wins over 1).

**Environment.** The child runs with `CUDA_DEVICE_ORDER=PCI_BUS_ID`, and so does
the `nvidia-smi` query, so `--gpu N` names the same board in both. `--gpu`
cannot be negative. A `CUDA_VISIBLE_DEVICES` of the form `GPU-<uuid>` is
accepted but cannot be checked against the "never GPU 0" rule, so the run says
so. `OMP_NUM_THREADS` and `MKL_NUM_THREADS` are recorded as inherited (a note
is printed when `OMP_NUM_THREADS` is unset).

The stage times are the `"timings"` of terco's `--json`
(`docs/sphinx/output.md`). The Fock number is the mean of builds 2..n; build 1
is its own row (`scf.fock_first`) because it pays one-off device setup. The
fock cases use `--no-incremental --thresh 1e-11`, as the cross-code harness
does, so the per-build time measures the integral engine and not convergence.

## The cases

`python3 run.py --list` prints them. Budget: at most 5 minutes per input on one
V100, one repeat. The sizes below were measured on a Tesla V100-DGXS-32GB
(2026-11-10); the time is the total wall time of one run.

| Case | Input | Guards | V100 |
|---|---|---|---|
| `fock_631g` | glycine chain x100, 6-31G, full builds | segmented SP shells, Cartesian | 105 s |
| `fock_631gs` | glycine chain x100, 6-31G* | SP shells plus d | 151 s |
| `fock_dz` | glycine chain x50, cc-pVDZ (quick) | s/p/d classes, screening | 70 s |
| `fock_dz_incr` | the same, default incremental Fock path | incremental schedule: Fock time per iteration and the iteration count | 56 s |
| `fock_tz` | water x40, cc-pVTZ, needs lmax3 | f classes | 148 s |
| `fock_tzvp` | glycine chain x20, def2-TZVP, needs lmax3 | f classes, mid-size cell | 176 s |
| `dft` | glycine chain x20, wB97X-V / def2-SVP | grid, XC, range-separated exchange, VV10 | 70 s |
| `df_w6_dz` | water hexamer (prism), cc-pVDZ `--df` (quick) | DF-JK | 2 s |
| `df_w6_tz` | water hexamer, cc-pVTZ `--df`, needs lmax3 | DF-JK with f | 3 s |
| `df_w6_atz` | water hexamer, aug-cc-pVTZ `--df`, needs lmax3 | DF-JK with diffuse functions | 3 s |
| `pcm` | glycine chain x20, def2-SVP IEF-PCM (quick) | PCM setup and per iteration | 19 s |
| `grad` | glycine chain x50, def2-SVP `--grad` (quick) | gradient kernels | 62 s |
| `rimp2` | glycine chain x20, cc-pVDZ / cc-pVDZ-RIFIT | RI-MP2 steps | 82 s |
| `dh` | glycine chain x20, B2PLYP / cc-pVDZ | double hybrid: four-centre SCF with XC, RI-MP2 for the PT2 part | 87 s |
| `ecp` | 24 CH3I, def2-SVP, iodine ECP | ECP matrix and SCF with an ECP | 46 s |

One repeat of everything takes about 18 minutes on one V100, so a full
comparison (a warm-up run and 3 repeats) takes a little over an hour and a
golden recording (a warm-up and 5 repeats) under two; the cases can be spread
over several GPUs with `--cases` and `--gpu`. `--quick` (`fock_dz`, `df_w6_dz`,
`pcm`, `grad`) takes about 10 minutes with the warm-up and 3 repeats.

Why these sizes. Density fitting pays below about a thousand basis functions,
so the DF cases are water hexamers and their stages are milliseconds; the
judging floor is lowered for them. Three cases are bounded by the memory of a
32 GB card and not by the time budget: PCM on the x50 glycine chain fails to
allocate (14 GB), the ECP caches grow as shells x ECP centres (32 CH3I already
fail, 64 would need 60 GB), and wB97X-V on the x50 chain takes 386 s, 283 s of
it VV10.

Geometries are in `systems/`, and every file there is used by a case. The
water cluster (`water40`) and glycine chains (`gly20`, `gly50`, `gly100`) are
copies of the inputs of the cross-code harness (`terco-fock-bench`) that the
Fock-build benchmarks use, so the bench stays comparable with those numbers.
`w6_prism.xyz` is the prism water hexamer of the maintainer's
`interesting_molecules` collection (MIT). `ch3i_24.xyz` is 24 copies of the
CH3I validated in `test/check_ecp_matrix.F90`, the first 24 sites of a 4x4x4
lattice, all parallel, with H...H contacts of 4.2 A and I...I contacts of
6.0 A. It is not a relaxed structure; a real molecule with one or two heavy
atoms and a large organic part would be a better ECP case.

## Recording a golden on a new architecture

1. Build terco for that machine (`-DTERCO_GPU_ARCH=...`; `-DTERCO_LMAX=3` if
   you want the f cases) and check the box is otherwise idle.
2. `python3 bench/regression/run.py --terco build/app/terco --gpu N --record --repeats 5`
   (5 or more: only then does the golden carry a usable spread). It writes
   `goldens/<arch>.csv` and prints what changed versus the previous golden,
   row by row; nothing is replaced silently. Cases that errored are not
   recorded (their old rows stay), and a partial `--record --cases a,b` keeps the
   other cases' rows and their own provenance.
3. Commit that CSV in a PR and say what machine, driver and build it came from.

The file is named by the GPU **family** slug (`v100.csv`, `a100.csv`,
`h100.csv`, `h200.csv`, `gh200.csv`; `cpu-<model>.csv` for a host build);
`--arch` overrides it, so a maintainer who wants one golden per host can use
`v100-<host>.csv`. Every row records the full `gpu_name` nvidia-smi reports
(for example `Tesla V100-SXM2-32GB`) and the host. At compare time a golden
whose `gpu_name` or `cpu_model` differs from this machine's draws a warning
(same family, different board or host: the judged device stages are still
compared, the host-side INFO stages will differ).

A golden is updated **only for an intended change**, stated in the PR: a
speed-up you made on purpose, a deliberate change of the numerics, a new
machine. Never to make a regression go away.

## The golden CSV

Flat, one row per `(case, stage)`, in the order of the cases in `cases.toml`
and then by stage name. Columns: `case, stage, median_s, min_s, max_s,
repeats, energy_hartree, iterations, gradient_rms`, then the provenance
`binary_commit, checkout_commit, date, host, cpu_model, gpu_name, driver, cuda,
omp_threads, mkl_threads, arch, build, compiler`. Times are fixed `%.6f`
seconds, the energy `%.10f`, the gradient RMS `%.10e`. The energy, iteration
count and gradient RMS are those of the case, repeated on its rows. For
`--rimp2` the energy is the SCF total plus the MP2 correlation.

`binary_commit` is the `git describe --always --dirty` that `terco --version`
reports for the **binary** (`commit=`), captured when CMake last configured, so
a later pull without a reconfigure leaves it stale; `unknown` outside a git
checkout. `checkout_commit` is the checkout `run.py` ran in. `driver` and
`cuda` are what `nvidia-smi` reports (the CUDA version is the driver's, not the
toolkit's).

## How the thresholds were chosen

`cases.toml` defaults: `fail` 1.10, `warn` 1.05, `min_seconds` 0.05,
`energy_tol` 1e-8. Once a golden carries 5 or more repeats the script derives
the fail ratio per stage as described above (`max(configured, 1 + 4 x spread)`).
Per-stage overrides go in a `[thresholds."<stage>"]` table and per-case ones in
the case (a case may also set `fail`, `warn` or `min_seconds` for all of its
stages).

Measured run-to-run spread of the judged stages on a V100 (Tesla
V100-DGXS-32GB, NVHPC 26.5, 5 repeats after one discarded warm-up run,
2026-11-10; the numbers of `goldens/v100.csv`):

| case | judged stage | median (s) | spread, (max - min)/median |
|---|---|---|---|
| `fock_631g` | `scf.fock_mean_after_first` | 6.088 | 0.05 % |
| `fock_631gs` | `scf.fock_mean_after_first` | 7.203 | 0.73 % |
| `fock_dz` | `scf.fock_mean_after_first` | 3.937 | 0.19 % |
| `fock_dz_incr` | `scf.fock_per_iteration` | 2.674 | 0.13 % |
| `fock_tz` | `scf.fock_mean_after_first` | 11.56 | 0.04 % |
| `fock_tzvp` | `scf.fock_mean_after_first` | 11.73 | 0.08 % |
| `dft` | `dft.grid` | 2.543 | 0.34 % |
| `dft` | `dft.vv10` | 45.36 | 0.25 % |
| `dft` | `dft.xc` | 3.421 | 5.06 % |
| `dft` | `scf.fock_mean_after_first` | 0.9632 | 0.35 % |
| `df_w6_dz` | `dfjk.step.k_accum_gemm` | 0.01501 | 0.48 % |
| `df_w6_dz` | `dfjk.step.k_half_transform` | 0.007002 | 1.94 % |
| `df_w6_dz` | `dfjk.step.mu_nu_p` | 0.05071 | 1.40 % |
| `df_w6_dz` | `scf.fock_mean_after_first` | 0.004946 | 2.28 % |
| `df_w6_tz` | `dfjk.step.k_accum_gemm` | 0.03063 | 0.19 % |
| `df_w6_tz` | `dfjk.step.k_half_transform` | 0.03157 | 0.36 % |
| `df_w6_tz` | `dfjk.step.mu_nu_p` | 0.1427 | 0.51 % |
| `df_w6_tz` | `scf.fock_mean_after_first` | 0.01336 | 0.18 % |
| `df_w6_atz` | `dfjk.step.k_accum_gemm` | 0.07759 | 0.16 % |
| `df_w6_atz` | `dfjk.step.k_half_transform` | 0.08392 | 0.37 % |
| `df_w6_atz` | `dfjk.step.mu_nu_p` | 0.2016 | 0.63 % |
| `df_w6_atz` | `scf.fock_mean_after_first` | 0.02991 | 0.13 % |
| `pcm` | `pcm.per_iteration` | 0.2746 | 0.10 % |
| `pcm` | `pcm.setup` | 2.855 | 0.79 % |
| `pcm` | `scf.fock_mean_after_first` | 0.5915 | 0.09 % |
| `grad` | `gradient.total` | 13.9 | 0.35 % |
| `rimp2` | `rimp2.step.cholesky` | 0.01836 | 4.32 % |
| `rimp2` | `rimp2.step.k_i_gemms` | 58.56 | 0.11 % |
| `rimp2` | `rimp2.step.mu_nu_p` | 1.584 | 0.38 % |
| `rimp2` | `rimp2.step.transform_ia_p` | 1.63 | 0.18 % |
| `dh` | `dft.xc` | 3.534 | 0.30 % |
| `dh` | `rimp2.step.k_i_gemms` | 58.55 | 0.11 % |
| `dh` | `rimp2.step.mu_nu_p` | 1.586 | 0.30 % |
| `dh` | `rimp2.step.transform_ia_p` | 1.63 | 0.34 % |
| `dh` | `scf.fock_mean_after_first` | 0.8301 | 0.12 % |
| `ecp` | `scf.fock_mean_after_first` | 0.5033 | 0.11 % |
| `ecp` | `scf.setup` | 36.22 | 0.18 % |

The Fock-build means repeat within 0.8 %, the RI-MP2 and DF steps within 2 %,
and the two noisy stages are the semilocal XC time of the wB97X-V case (5 %)
and the 18 ms RI-MP2 Cholesky step (4 %). So the configured 1.10 applies to
almost every stage, and the spread rule widens it only for those two (to about
1.20). Rerunning four cases alone on one GPU reproduced the golden, which was
recorded with three GPUs busy at once, within 0.3 %. A job contending for the
same GPU is reported as a regression (PCM case: 1.8x to 2.6x on its judged
stages, energy unchanged). Fill in the same table when a golden is recorded on
another architecture.

## Checking the checker

`python3 bench/regression/run.py --self-test` needs no terco: it tests the
parsing, the median and verdict logic (at the warn and fail ratios, the spread
threshold, SMALL growth), the golden round trip and row order, the GPU-0
refusal and, with a fake `terco`, the whole record / compare / failure /
timeout path: per-repeat energy failures, errored cases not aborting the
others, warm-up runs discarded, INFO stages never failing, build-mismatch
refusal, a partial record, a relative `--out`, and the timeout killing the
process group.
