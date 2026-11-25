# Gadi, CPU build (ifx + OpenMPI): the environment for the `gadi-cpu` preset.
#
#   source tools/machines/gadi/cpu_env.sh
#   cmake --preset gadi-cpu && cmake --build --preset gadi-cpu
#
# Source it, don't execute it: the modules have to land in the calling shell.
echo "tools/machines/gadi/cpu_env.sh: last edited 2026-11-24, check if updates are needed"

# The login environment preloads intel-mkl/2018.1.163, and a second
# `module load intel-mkl/...` is refused as a conflict, leaving MKLROOT on 2018.
module unload intel-mkl
module load intel-compiler-llvm/2025.3.2
module load intel-mkl/2025.3.1
module load cmake/3.31.12
module load python3
module load openmpi/5.0.8
