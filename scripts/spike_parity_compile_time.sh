#!/bin/sh
# PARITY-SKIP SPIKE: host compile time of the masked kernel variants.
# usage: spike_parity_compile_time.sh <build dir with modules/> <python> <scratch dir> [reps]
# Emits each (class, mask) kernel as its own module with
# scripts/spike_parity_gen.py and times `gfortran-15 -O3 -c` on it (best of
# `reps`, default 3) of CPU seconds (user+sys, so a busy machine matters
# less), the flags of the Release build. Nothing is written
# outside the scratch dir.
set -e
BUILD=$1; PY=$2; OUT=$3; REPS=${4:-3}
HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"
printf '%-8s %-5s %8s %10s %10s\n' class mask lines seconds objbytes
for spec in 1111:0 1111:4 1111:6 2020:0 2020:4 2121:0 2121:4 2121:6 2211:0 2211:4 2222:0 2222:4 2222:6; do
  tag=${spec%%:*}; mask=${spec##*:}
  f="$OUT/k${tag}z${mask}.F90"
  "$PY" "$HERE/spike_parity_gen.py" module "$spec" --name "spk_${tag}_${mask}" -o "$f"
  best=$("$PY" - "$f" "$OUT" "$BUILD" "$REPS" <<'PYEOF'
import resource, subprocess, sys
f, out, build, reps = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
best = 1e9
for _ in range(reps):
    r0 = resource.getrusage(resource.RUSAGE_CHILDREN)
    subprocess.run(["/opt/homebrew/bin/gfortran-15", "-O3", "-J", out, "-I", build + "/modules", "-I", build + "/_deps/pic-build/modules", "-c", f, "-o", out + "/k.o"], check=True)
    r1 = resource.getrusage(resource.RUSAGE_CHILDREN)
    best = min(best, r1.ru_utime - r0.ru_utime + r1.ru_stime - r0.ru_stime)
print(round(best, 2))
PYEOF
)
  printf '%-8s %-5s %8s %10s %10s\n' "$tag" "$mask" "$(wc -l < "$f")" "$best" "$(stat -f %z "$OUT/k.o")"
done
