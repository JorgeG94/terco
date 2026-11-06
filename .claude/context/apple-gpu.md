# Apple GPU feasibility (S0b)

Date of research: 2026-11-06. Research only; no product code or CMake changed.

## Question

Can terco's `do concurrent` kernels (all `real(dp)`, double precision) run on
the M5 Pro GPU while staying standard Fortran?

## Conclusion: NO-GO

Two independent blockers, either sufficient:

1. Metal has no fp64. terco's kernels are `real(dp)`; they cannot be
   expressed in MSL, and the Metal feature tables list no fp64 capability
   for any Apple GPU family.
2. No standard-Fortran toolchain targets Apple GPUs. flang's
   `do concurrent` offload goes through OpenMP target offload, whose LLVM
   plugins are NVPTX and AMDGPU (plus experimental SPIR-V), not Metal.

Revisit only if Apple ships fp64 in MSL or the project accepts fp32 (which
breaks the 1e-8-level accuracy gates), or an emulated-double approach is
accepted as a separate, non-standard-Fortran code path.

## Findings

### F1. MSL does not support `double`
Metal Shading Language Specification v4.1, section 2.1 "Scalar Data Types":
"Metal does not support the double, long long, unsigned long long, and long
double data types." (Read from the PDF fetched 2026-11-06, page 24 of 383.)
https://developer.apple.com/metal/Metal-Shading-Language-Specification.pdf

### F2. Feature-set tables list 64-bit integer, not fp64
Metal Feature Set Tables list "64-bit integer math" (Apple3 and later) and
64-bit atomics (Apple9), with no floating-point 64 row. Searching the PDF
text for "double", "float64" and "fp64" finds nothing.
https://developer.apple.com/metal/Metal-Feature-Set-Tables.pdf

### F3. Corroboration (secondary)
PyTorch MPS rejects float64 tensors; Apple forum threads and an AdaptiveCpp
discussion report double failing in the shader compiler; the only
workaround is integer-based emulation, which is slow.
https://developer.apple.com/forums/thread/709663
https://github.com/orgs/alicevision/discussions/2160

### F4. flang `-fdo-concurrent-to-openmp=device`
Maps `do concurrent` to `omp target teams distribute parallel do`; needs
`-fopenmp` plus offload target flags. Status (LLVM 22 docs): "basic support
that is still in progress", tested on 1D/2D saxpy; no locality specifiers or
reductions yet. No target list is given.
https://flang.llvm.org/docs/DoConcurrentConversionToOpenMP.html

### F5. LLVM OpenMP offload targets
The OpenMP offload docs describe plugins for NVPTX (nvptx64-nvidia-cuda)
and AMDGPU (amdgcn-amd-amdhsa). SPIR-V and Apple/Metal are not mentioned.
https://openmp.llvm.org/SupportAndFAQ.html

### F6. Hands-on test with Homebrew flang 23.1.5 (arm64-apple-darwin25.5.0)
`brew install flang` worked (no sudo). Scratch program
`/private/tmp/claude-501/s0b-scratch/dc.f90` (1e6-element `do concurrent`,
real(8) a = 3b+1, prints 7.). Results:
- plain `flang -O2`: compiles, runs, prints 7.
- `-fopenmp -fdo-concurrent-to-openmp=host`: compiles, runs, prints 7.
- `-fopenmp -fdo-concurrent-to-openmp=device` with no offload target:
  compiles silently (no device code generated; runs on host).
- with `--offload-arch=gfx90a`: fails, "cannot find ROCm device library".
- with `--offload-targets=spirv64`: fails, no `libomptarget-spirv.bc` in
  the Homebrew LLVM lib directory.
- the Homebrew LLVM ships only `libomp.dylib`; no device runtime plugins.
So the host path works on this Mac; there is no way to name a Metal target.

### F7. Other routes
- SPIR-V to Metal (SPIRV-Cross, MoltenVK): the output is still MSL, so fp64
  is still missing (F1). MoltenVK documentation lists no 64-bit float
  extension.
  https://github.com/KhronosGroup/MoltenVK/blob/main/Docs/MoltenVK_Runtime_UserGuide.md
- Hand-written Metal kernels, called from Fortran via C interop: possible
  in fp32 only, and it abandons the standard-Fortran premise (separate
  kernel language, Objective-C/Swift host glue, duplicated kernels).
- Emulated double (double-float / integer based): no standard Fortran path,
  large slowdown, and it needs the same hand-written Metal route.

## Not run
No GPU code ran. Build and test definition-of-done items do not apply to
this research stage.
