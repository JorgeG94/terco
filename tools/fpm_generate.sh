#!/usr/bin/env bash
#
# Generate the kernels fpm cannot generate itself, then let fpm build.
#
# The CMake build generates the Rys kernels (energy and gradient, one module
# per class, plus their dispatchers) into the build tree at configure time.
# fpm cannot run a generator, so this script does the same generation into
# src/generated/fpm/, which fpm then compiles along with everything else
# under src/. Run it once before `fpm build`, and again after changing
# scripts/.
#
#   tools/fpm_generate.sh
#   fpm build --compiler mpifort
#
# The committed kernels in src/generated are TRC_LMAX=2, lexicographic
# Cartesian order, so this generates for the same configuration. Any other
# ceiling or order needs the CMake build. src/generated/fpm/ is gitignored and
# never seen by CMake (its glob of src/generated is not recursive).
#
# Environment:
#   TERCO_RYS_UNROLL_MAX  as the CMake option of the same name (default 81)
#   JOBS                  parallel generator processes (default: all cores)
#
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
lmax=2
unroll="${TERCO_RYS_UNROLL_MAX:-81}"
jobs="${JOBS:-$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4)}"
out="${root}/src/generated/fpm"
py="${PYTHON:-python3}"

tmp="$(mktemp -d)"
trap 'rm -rf "${tmp}"' EXIT
stage="${tmp}/stage"
template="${tmp}/hgp_template"
mkdir -p "${stage}" "${template}"

# Every class up to lmax; the energy Rys set is every class with a d shell
# plus (pp|pp) and (pp|ps), exactly as src/CMakeLists.txt selects it. The
# gradient set is every class up to d.
all=()
rys=()
for la in $(seq 0 ${lmax}); do
  for lb in $(seq 0 ${lmax}); do
    for lc in $(seq 0 ${lmax}); do
      for ld in $(seq 0 ${lmax}); do
        tag="${la}${lb}${lc}${ld}"
        all+=("${tag}")
        if [[ "${tag}" =~ [2-9] || "${tag}" == 1111 || "${tag}" == 1110 ]]; then
          rys+=("${tag}")
        fi
      done
    done
  done
done
join() { local IFS=,; echo "$*"; }

cd "${root}"
echo "fpm_generate: HGP template for ${#all[@]} classes"
"${py}" scripts/gen_perclass.py --lmax ${lmax} --classes "$(join "${all[@]}")" \
  --skeleton --split "${template}" > /dev/null

echo "fpm_generate: ${#rys[@]} Rys energy kernels, ${#all[@]} gradient kernels (${jobs} jobs)"
{
  for tag in "${rys[@]}"; do
    echo "scripts/gen_rys.py ${tag} --src-dir ${template} --unroll-max-comp ${unroll} -o ${stage}/trc_pc_rys${tag}.F90"
  done
  for tag in "${all[@]}"; do
    echo "scripts/gen_rys_grad.py ${tag} --src-dir ${template} --unroll-max-comp ${unroll} -o ${stage}/trc_pc_rysg${tag}.F90"
  done
} | xargs -P "${jobs}" -I{} sh -c "\"${py}\" {} > /dev/null"

"${py}" scripts/gen_rys_dispatch.py --classes "$(join "${rys[@]}")" \
  -o "${stage}/trc_rys_kernels.F90" > /dev/null
"${py}" scripts/gen_rys_grad_dispatch.py --classes "$(join "${all[@]}")" \
  -o "${stage}/trc_rys_grad_kernels.F90" > /dev/null

# Swap in only once everything generated, so a failure never leaves a
# half-written set behind for fpm to pick up.
rm -rf "${out}"
mv "${stage}" "${out}"
echo "fpm_generate: $(ls "${out}" | wc -l) files in src/generated/fpm"
