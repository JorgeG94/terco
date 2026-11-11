# Goldens

One CSV per architecture, named by the GPU family slug `run.py` derives from
`nvidia-smi` (overridable with `--arch`, for example `v100-<host>`) or, for a
host build, from the CPU model. The rows record the full GPU name and the host:

| File | Machine |
|---|---|
| `v100.csv` | Tesla V100 |
| `a100.csv` | A100 |
| `h100.csv` | H100 |
| `h200.csv` | H200 |
| `gh200.csv` | GH200 (Grace Hopper) |
| `cpu-<model>.csv` | a host build, e.g. `cpu-intel-xeon-gold-6248.csv` |

None is committed yet: they hold real timings, which are recorded on the
machine they describe (`run.py --record`, see ../README.md) and added by a PR
that says why. Do not write numbers into these files by hand.
