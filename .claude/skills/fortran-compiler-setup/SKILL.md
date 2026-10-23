---
name: fortran-compiler-setup
description: Install recent Fortran compilers (LLVM flang, gfortran 13-16, LFortran, Intel ifx) in a fresh Claude cloud sandbox (Ubuntu 24.04 behind an egress allowlist), including which network domains must be enabled. Use before compiling or testing Fortran code with any of these compilers, or when a compiler download fails with a proxy 403.
---

# Fortran compilers in the Claude cloud sandbox

Each cloud session starts from a fresh Ubuntu 24.04 x86-64 container (2 vCPUs sharing one core),
so compilers are reinstalled every session. All outbound HTTPS goes through an allowlisting proxy.

## Ground rules

1. **Binary first.** Use a prebuilt binary whenever one is reachable; build from source only when no
   binary route works. Cost of a source build: LLVM flang — never do it here; GCC — only on explicit
   request (coarrays, offloading); LFortran — affordable (15–25 min).
2. **Never patch or rebuild an installed compiler** to diagnose or work around a bug. The goal is a
   minimal reproducer reported to the vendor, not a fix.
3. Start long installs in the background first thing, then do other work.
4. Write PATH settings to `~/fortran-env.sh` and `source ~/fortran-env.sh` in each command:
   `HOME=/root` and `~/.bashrc` exits early for non-interactive shells, so `.bashrc` edits don't take effect.
   The scripts below append to that file.

## Step 0 — check the network

```sh
bash scripts/check_network.sh
```
Prints OK/BLOCKED per host and the proxy's recent denials. A proxy 403 is a policy denial: don't retry or
route around it — the domain has to be added to the environment's network allowlist by the user.

| Needed for | Domains to allow |
|---|---|
| flang, xPack gfortran, LFortran release tarball, pixi binary | `github.com`, `release-assets.githubusercontent.com` (asset CDN; older links use `objects.githubusercontent.com`) |
| LFortran from git main | `github.com` (git over HTTPS) |
| apt packages (gfortran 13/14, LLVM 20 dev libs, OpenCoarrays) | `archive.ubuntu.com`, `security.ubuntu.com` |
| GCC 16 source build (coarrays) | `gcc.gnu.org` (release tarball and `contrib/download_prerequisites`) |
| Pixi / conda-forge binaries (LFortran, gfortran 15.3/16.2, flang) | `conda.anaconda.org`, `prefix.dev`, `repo.prefix.dev`, `fast.prefix.dev`; `pixi.sh` only for the installer script (not needed, see below) |
| Intel ifx | `apt.repos.intel.com`, `registrationcenter-download.intel.com` |
| Optional | `raw.githubusercontent.com` (single files), `pypi.org`/`files.pythonhosted.org` (newer CMake via pip; normally allowed by default) |

GitHub *web pages* (`/releases/latest`), `api.github.com` and `codeload.github.com` are usually blocked
even when release downloads work. To find the newest version, list tags with git instead:
```sh
git ls-remote --tags https://github.com/llvm/llvm-project 'llvmorg-2*' | grep -v -e '\^{}' -e rc | awk -F/ '{print $3}' | sort -V | tail -3
git ls-remote --tags https://github.com/xpack-dev-tools/gcc-xpack | awk -F/ '{print $3}' | grep -v '\^{}' | sort -V | tail -3
git ls-remote --tags https://github.com/lfortran/lfortran 'v0.*' | awk -F/ '{print $3}' | grep -v '\^{}' | sort -V | tail -3
```

## Which route for which compiler

| Compiler | Route | Time | Script |
|---|---|---|---|
| **flang 23.1.3** | LLVM GitHub release tarball (binary) | ~3 min, ~3 GB on disk | `install_flang.sh [23.1.3]` |
| **gfortran 13.3 / 14.2** | Ubuntu apt (binary) | ~1 min | `install_gfortran.sh apt` |
| **gfortran 15.2** | xPack GitHub release (binary) | <1 min | `install_gfortran.sh xpack [15.2.0-1]` |
| gfortran 15.3 / 16.2 | conda-forge via pixi (binary; needs conda domains; **untested** here) | ~2 min | `install_gfortran.sh conda 16.2` |
| **gfortran 16.2 + caf_shmem coarrays** | GCC source build — **on demand only** | ~30–75 min | `install_gfortran.sh source [16.2.0]` |
| **LFortran 0.67.0** | conda-forge binary via pixi (needs conda domains; **untested** here) | ~2 min | `install_lfortran.sh pixi` |
| **LFortran git main** | source, apt LLVM 20, Debug + stack traces | ~15–25 min | `install_lfortran.sh main` |
| LFortran release | source tarball, apt LLVM 20, Release | ~25 min | `install_lfortran.sh release [0.67.0]` |
| LFortran git main, LLVM 11 | upstream `pixi run build` (CI reference LLVM; needs conda domains) | ~20 min | `install_lfortran.sh pixi-main` |
| **Intel ifx 2026.2** | Intel oneAPI apt repo (binary), **>2 GB** | ~3 min | `install_ifx.sh [2026.2]` |

Versions current as of 2026-10-19; probe for newer tags before installing.

## flang (LLVM)

Use the upstream release tarball — the right approach. Don't use Ubuntu's `flang-18` (too old: no `class(*)`
polymorphism). The script extracts only flang, lld and the libraries from the ~1 GB tarball:
```sh
mkdir -p ~/llvm23 && cd ~/llvm23
curl -fsSL https://github.com/llvm/llvm-project/releases/download/llvmorg-23.1.3/LLVM-23.1.3-Linux-X64.tar.xz \
  | tar xJ --wildcards '*/bin/flang*' '*/bin/ld.lld' '*/bin/lld' '*/lib/*' '*/include/flang*'
```
- `libomp.so` is in `lib/x86_64-unknown-linux-gnu/` (the script adds it to `LD_LIBRARY_PATH`).
- flang's default OpenMP version is 3.1; pass `-fopenmp-version=52` for 5.x clauses.
- `-flto` needs `ld.lld` from the same tarball (no LLVMgold.so is shipped).
- conda-forge also has `flang 23.1.3`; the GitHub tarball needs fewer domains.

## gfortran

- **apt:** `gfortran` (13.3) and `gfortran-14` (14.2). Run `apt-get update` first; the image's lists are stale.
- **xPack 15.2:** link with `-Wl,-rpath,$GCC_XPACK/lib64` so the matching `libgfortran` is found at run time.
  xPack has no GCC 16 build (as of Oct 2026); the ubuntu-toolchain PPA is blocked.
- **conda-forge 15.3 / 16.2** (`pixi global install gfortran=16.2`): saves a source build if the conda domains
  are enabled. Untested in the sandbox; check whether it ships `libcaf_shmem` before relying on it for coarrays.

### Optional: gfortran 16 with shared-memory coarrays (caf_shmem) — build on demand only

GCC 16 adds `libcaf_shmem`, a multi-image coarray runtime using forked processes on one node (no MPI).
No reachable binary includes it, so it needs a genuine GCC build. **Only start it when the task needs
coarrays**, and start it immediately in the background. Expected time: roughly 30 min on a core with two
threads by the user's estimate; the logged October 2026 sandbox build took ~75 min, so budget for that.

```sh
bash scripts/install_gfortran.sh source 16.2.0      # configure, then make in background
tail -f ~/gcc-src/build.log                          # done when it prints: GCC BUILD DONE
```
It configures with `--enable-languages=c,c++,fortran --disable-bootstrap --disable-multilib --disable-libsanitizer`.
Usage once built:
```sh
source ~/fortran-env.sh
gfortran -fcoarray=lib scripts/caf_smoke.f90 -lcaf_shmem -o caf_smoke    # or -lcaf_single for one image
GFORTRAN_NUM_IMAGES=4 ./caf_smoke                                        # prints: OK: 4 images
```
- Every object in a program must use the same `-fcoarray` mode.
- Images are fork()ed processes: remote access to a pointer component (`z[p]%ptr`) works only when the target is
  a coarray or an allocatable component of one; private heap targets segfault.
- caf_shmem in GCC 16 has known errors that can make processes hang; fixes are expected in GCC 17. Use a
  timeout (`timeout 60 ./a.out`) and treat hangs as possible runtime bugs, not necessarily user errors.
- A second runtime for cross-checking: OpenCoarrays 2.9.2 from apt (`libcoarrays-openmpi-dev`, gfortran 13), with
  `OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 OMPI_MCA_osc=pt2pt`, then `caf` / `cafrun -n 4 --oversubscribe`.
- The 2 vCPUs share one core's CPU time: correctness testing is fine, parallel scaling can't be measured here.

## LFortran

**Pixi works only if the conda-forge domains are enabled** (`conda.anaconda.org`, `prefix.dev`, `repo.prefix.dev`,
`fast.prefix.dev`). The pixi binary itself installs from GitHub releases (`install_pixi.sh`) — the `pixi.sh`
installer script is not needed. With those domains:
- `install_lfortran.sh pixi` → conda-forge **binary** of the latest release (0.67.0, built against LLVM 21).
  This is the binary-first route; it needs a system C compiler for linking (`build-essential`).
- `install_lfortran.sh pixi-main` → upstream's own workflow (`pixi run build`, default env LLVM 11 =
  LFortran's CI reference), compiler at `build/llvm11/src/bin/lfortran`.

Without them (the default allowlist), build from source against apt LLVM 20 (`install_lfortran.sh main`):
```sh
apt-get install -y llvm-20-dev clang-20 libclang-rt-20-dev re2c bison zlib1g-dev libzstd-dev libunwind-dev \
                   binutils-dev cmake ninja-build python3 python-is-python3
git clone https://github.com/lfortran/lfortran.git && cd lfortran
export CMAKE_PREFIX_PATH_LFORTRAN=/usr/lib/llvm-20 PATH=/usr/lib/llvm-20/bin:$PATH
bash ./build1.sh build -DLLVM_DIR=/usr/lib/llvm-20/lib/cmake/llvm > build.log 2>&1   # binary: build/src/bin/lfortran
```
- `libunwind-dev` is required (configure fails without it). `build1.sh` runs the re2c/bison code generation and
  needs full git history for `git describe` — use a full clone, not `--depth 1`.
- **Bug reports should be made against git main** — releases lag main considerably. Use the binary for quick
  triage, then confirm on main before filing. Before filing, re-check on LLVM 11 (`pixi run -e llvm11`) if
  available, since upstream reference outputs are defined against it.
- CMake ≥ 3.31 knows the LFortran compiler ID (Ubuntu ships 3.28): `python3 -m venv venv && venv/bin/pip install cmake`.
- Legacy code flags: `--fixed-form`, `--implicit-typing`, `--implicit-interface`, `--legacy-array-sections`,
  `--separate-compilation` (must be used for all files). `--std=legacy` forces fixed-form parsing even for `.f90`
  files; use `--std=f23` for free-form files in mixed projects.

## Intel ifx

Requires `apt.repos.intel.com` and `registrationcenter-download.intel.com`. **Large: >2 GB** (about 0.85 GB of
packages, 2.4 GB installed for 2026.2); ~3 min.
```sh
curl -fsSL https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
  | gpg --dearmor --yes -o /usr/share/keyrings/oneapi-archive-keyring.gpg
echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" \
  > /etc/apt/sources.list.d/oneAPI.list
apt-get update && apt-get install -y intel-fortran-compiler-2026.2   # apt-cache search '^intel-fortran-compiler-20'
source /opt/intel/oneapi/setvars.sh                                   # ifx 2026.2.0
```
- `apt-get update` may warn about other blocked repos (e.g. docker); that's harmless.
- CMake compiler ID is `IntelLLVM`; OpenMP flag `-qopenmp` / `-fiopenmp`. PyPI only has the runtimes.

## Before filing a compiler bug

- Reduce to a minimal standalone reproducer; record the exact `--version` output and flags.
- Cross-check against at least one other compiler to separate compiler bugs from non-conforming code.
- Don't read the compiler's source to find the cause and don't patch it; describe the observed vs. expected
  behaviour and move on.
