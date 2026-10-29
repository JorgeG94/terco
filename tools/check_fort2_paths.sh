#!/bin/bash
#
# nvfortran's fort2 has segfaulted (signal 11, or "ili_hash_func: invalid opc
# 65535") compiling generated kernels for the GPU, on classes with two or
# more s shells, at SOME absolute build-path lengths and not others: whether
# it fires depends on fort2's heap layout, which the path length moves (see
# _jonly_digestion in scripts/gen_perclass.py). A build on one machine says
# nothing about another; Gadi hit it on rys0201/0210/1021 where the V100 box
# did not. So compile the at-risk objects for the GPU at a sweep of path
# lengths. No GPU is needed: this only compiles.
#
#   tools/check_fort2_paths.sh [arch] [min_len] [max_len] [step]
#     arch: TERCO_GPU_ARCH (default cc90); lengths of the build directory's
#     absolute path (default 45..200 step 15). nvfortran must be on PATH.
#
# Exit status 1 if any object fails at any length; the failing (length,
# object) pairs are listed.
set -u
SRC=$(cd "$(dirname "$0")/.." && pwd)
ARCH=${1:-cc90}
LMIN=${2:-45}; LMAX=${3:-200}; STEP=${4:-15}
ROOT=${RUNNER_TEMP:-${TMPDIR:-/tmp}}/fort2paths
mkdir -p "$ROOT"
JOBS=${JOBS:-$(nproc)}
fails=()
for len in $(seq "$LMIN" "$STEP" "$LMAX"); do
  pad=$(( len - ${#ROOT} - 1 ))
  if [ $pad -lt 1 ]; then echo "skip length $len (root $ROOT is ${#ROOT} long)"; continue; fi
  dir="$ROOT/$(printf 'b%.0s' $(seq 1 $pad))"
  rm -rf "$dir"
  if ! cmake -G Ninja -S "$SRC" -B "$dir" -DCMAKE_BUILD_TYPE=Release \
       -DCMAKE_Fortran_COMPILER=nvfortran -DTERCO_GPU_ARCH="$ARCH" \
       -DTERCO_ENABLE_TESTING=OFF -DTERCO_BUILD_APP=OFF > "$dir.cfg.log" 2>&1; then
    echo "length $len: configure failed (see $dir.cfg.log)"; fails+=("$len:configure"); continue
  fi
  # At-risk objects: HGP (pc_k) and Rys class kernels with two or more s shells.
  objs=$(ninja -C "$dir" -t targets all 2>/dev/null \
         | grep -oE 'CMakeFiles/terco\.dir/[^:]*(trc_pc_k|trc_pc_rys)[0-9]{4}\.F90\.o' \
         | sort -u | while read -r o; do
             tag=$(echo "$o" | grep -oE '[0-9]{4}\.F90\.o' | cut -c1-4)
             nz=$(echo "$tag" | tr -cd 0 | wc -c)
             [ "$nz" -ge 2 ] && echo "$o"
           done)
  n=$(echo "$objs" | grep -c . || true)
  if ninja -C "$dir" -k 0 -j "$JOBS" $objs > "$dir.build.log" 2>&1; then
    echo "length $len: $n objects OK"
  else
    bad=$(grep -oE 'FAILED: [^ ]*\.o' "$dir.build.log" | sed 's/FAILED: //' | xargs -n1 basename | tr '\n' ' ')
    echo "length $len: FAILED $bad"
    grep -m3 -E 'fort2|signal|ili_hash' "$dir.build.log" | sed 's/^/    /'
    fails+=("$len:$bad")
  fi
done
if [ ${#fails[@]} -gt 0 ]; then
  echo "fort2 path-length sweep FAILED:"; printf '  %s\n' "${fails[@]}"; exit 1
fi
echo "fort2 path-length sweep passed ($LMIN..$LMAX step $STEP, arch $ARCH)"
