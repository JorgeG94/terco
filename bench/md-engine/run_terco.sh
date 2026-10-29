#!/bin/bash
# terco per-Fock timings for the MD engine.
#   B=<build dir> GPU=<device> ./run_terco.sh [nsys]
# Runs w5 and w20 cc-pVQZ RHF for 4 iterations (iteration 4 is the timed
# full build), with the 17-digit pyscf cc-pVQZ in this directory, and
# prints energy and fock(s) of iteration 4. With "nsys", also writes
# <xyz>.sqlite for join_profile.py (needs nsys on PATH).
set -u
H=$(cd "$(dirname "$0")" && pwd)
B=${B:?set B to the build directory}
GPU=${GPU:-0}
for x in w5 w20; do
    log=$H/out_$x.log
    if [ "${1:-}" = nsys ]; then
      CUDA_VISIBLE_DEVICES=$GPU nsys profile -t cuda,nvtx -o $H/$x --force-overwrite true \
        $B/app/terco $H/$x.xyz -b $H/cc-pvqz-pyscf17.json --iterations 4 > $log 2>&1
      nsys export -t sqlite -o $H/$x.sqlite --force-overwrite true $H/$x.nsys-rep > /dev/null 2>&1
    else
      CUDA_VISIBLE_DEVICES=$GPU $B/app/terco $H/$x.xyz -b $H/cc-pvqz-pyscf17.json \
        --iterations 4 > $log 2>&1
    fi
    echo "$x: $(grep -E '^ +4 ' $log | awk '{print "E",$2," fock",$6}')"
done
