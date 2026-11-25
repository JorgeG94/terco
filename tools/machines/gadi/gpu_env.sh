# Gadi, GPU build (nvfortran + OpenMPI): the environment for the
# `gadi-gpu-cc70` (V100) and `gadi-gpu-cc90` (Hopper) presets.
#
#   source tools/machines/gadi/gpu_env.sh
#   cmake --preset gadi-gpu-cc70 && cmake --build --preset gadi-gpu-cc70
#
# Source it, don't execute it: the modules have to land in the calling shell.
echo "tools/machines/gadi/gpu_env.sh: last edited 2026-11-24, check if updates are needed"

# The login environment preloads intel-mkl/2018.1.163, and a second
# `module load intel-mkl/...` is refused as a conflict, leaving MKLROOT on 2018.
module unload intel-mkl
module load nvidia-hpc-sdk/26.9
module load openmpi/5.0.5
module load intel-mkl/2025.3.1
module load cmake/3.31.12
module load python3
module load cuda
