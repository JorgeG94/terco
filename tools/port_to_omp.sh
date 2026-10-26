#!/bin/sh
# Generate the OpenMP-offload tree (what portability/amd-omp used to carry by
# hand) from this checkout: nothing in the result is edited, so it can be
# thrown away and regenerated at any time.
#
#   tools/port_to_omp.sh DEST [gpu|cpu]
#
# DEST must not exist. Copies every tracked file of the working tree, then
# runs the two converters in their required order (see docs/GENERATORS.md)
# on src, test and app. src/xc/vendor is put back as upstream: it already
# carries both annotations, and the converter only strips the OpenACC line.
#
# Build the result with the frontier-gpu preset (TERCO_ENABLE_HIP); at a
# TERCO_LMAX other than the shipped one CMake converts the regenerated
# kernels itself.
set -eu
dest=${1:?usage: tools/port_to_omp.sh DEST [gpu|cpu]}
target=${2:-gpu}
here=$(cd "$(dirname "$0")/.." && pwd)
py=${PYTHON:-python3}
[ ! -e "$dest" ] || { echo "port_to_omp: $dest already exists" >&2; exit 1; }
mkdir -p "$dest"
(cd "$here" && git ls-files -z | xargs -0 cp --parents -t "$dest")
cd "$dest"
"$py" tools/dc_to_omp.py --target "$target" --write -j 0 src test app
"$py" tools/acc_to_omp.py --target "$target" --translate-atomics --write -j 0 src test app
(cd "$here" && cp -p src/xc/vendor/* "$dest/src/xc/vendor/")
