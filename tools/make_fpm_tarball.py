#!/usr/bin/env python3
"""Package a release for fpm at one TERCO_LMAX: terco-<version>-lmax<N>.tar.gz.

fpm cannot run terco's kernel generators, and the committed src/generated is
the TERCO_LMAX=2 set only, so an fpm user at any other ceiling had to run the
generators by hand. This does what a CMake configure does, once, and packs the
result: the tree at <ref> with src/generated replaced by exactly the generated
sources CMake compiles for that ceiling, and fpm.toml's macros set to CMake's
definitions for it. The receiver needs only fpm and a compiler:

    tar xf terco-0.1.0-lmax3.tar.gz && cd terco-0.1.0-lmax3
    fpm install --compiler mpifort --prefix ... --flag "..."

What CMake compiles is read from its file API (codemodel-v2), not from a copy
of src/CMakeLists.txt's class selection, so the two cannot drift. Generation
does not depend on the compiler or on TERCO_ENABLE_GPU, so a CPU configure
(with any Fortran compiler CMake accepts; nothing is compiled) serves a GPU
build.

    tools/make_fpm_tarball.py --lmax 3 [--ref v0.1.0] [--out dist]
"""
import argparse
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile

# Not shipped in the tarball: the tests' inputs and references (44 MB), the
# benchmark inputs, CI and agent notes. None of them is needed to build.
DROP = ["test", "bench", ".claude", ".github", "CLAUDE.md", "AGENTS.md"]
# Data files the CLI finds through the same lookup as a basis (the minao and
# SAP guesses, Hirshfeld's proatoms). They live with the tests in the repo.
DATA = ["minao.json", "sap_helfem_large.json", "hirshfeld_proatoms.json"]


def run(cmd, **kw):
    print("+", " ".join(cmd), flush=True)
    subprocess.run(cmd, check=True, **kw)


def codemodel_target(build, name):
    reply = os.path.join(build, ".cmake", "api", "v1", "reply")
    cm = json.load(open(glob.glob(os.path.join(reply, "codemodel-v2-*.json"))[0]))
    for t in cm["configurations"][0]["targets"]:
        if t["name"] == name:
            return json.load(open(os.path.join(reply, t["jsonFile"])))
    sys.exit(f"make_fpm_tarball: no target {name} in the CMake codemodel")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lmax", type=int, required=True)
    ap.add_argument("--ref", default="HEAD", help="git ref to package (the release tag)")
    ap.add_argument("--out", default="dist")
    ap.add_argument("--fortran", default=os.environ.get("FC", "gfortran"),
                    help="any Fortran compiler CMake accepts; nothing is compiled")
    ap.add_argument("--keep", action="store_true", help="keep the work directory")
    a = ap.parse_args()

    root = subprocess.run(["git", "rev-parse", "--show-toplevel"], check=True,
                          capture_output=True, text=True).stdout.strip()
    manifest = subprocess.run(["git", "-C", root, "show", f"{a.ref}:fpm.toml"], check=True,
                              capture_output=True, text=True).stdout
    version = re.search(r'^version\s*=\s*"([^"]+)"', manifest, re.M).group(1)
    name = f"terco-{version}-lmax{a.lmax}"
    out = os.path.abspath(a.out)
    os.makedirs(out, exist_ok=True)

    work = tempfile.mkdtemp(prefix="terco-fpm-")
    tree = os.path.join(work, name)
    build = os.path.join(work, "cmake")
    try:
        # The tree at the ref, generators included: generate from the same
        # scripts the tarball ships.
        archive = subprocess.run(["git", "-C", root, "archive", "--format=tar", f"--prefix={name}/", a.ref],
                                 check=True, capture_output=True).stdout
        tarfile.open(fileobj=__import__("io").BytesIO(archive)).extractall(work, filter="data")

        os.makedirs(os.path.join(build, ".cmake", "api", "v1", "query"))
        open(os.path.join(build, ".cmake", "api", "v1", "query", "codemodel-v2"), "w").close()
        run(["cmake", "-G", "Ninja", "-S", tree, "-B", build, "-DCMAKE_BUILD_TYPE=Release",
             f"-DCMAKE_Fortran_COMPILER={a.fortran}", "-DTERCO_ENABLE_GPU=OFF",
             "-DTERCO_ENABLE_TESTING=OFF", f"-DTERCO_LMAX={a.lmax}"],
            stdout=subprocess.DEVNULL)

        lib = codemodel_target(build, "terco")
        srcs = set()
        for s in lib["sources"]:
            p = s["path"] if os.path.isabs(s["path"]) else os.path.join(tree, s["path"])
            srcs.add(os.path.realpath(p))
        group = [g for g in lib["compileGroups"] if g["language"] == "Fortran"]
        if len(group) != 1:
            sys.exit("make_fpm_tarball: expected one Fortran compile group (a per-file definition fpm "
                     "cannot express?)")
        defines = [d["define"] for d in group[0].get("defines", []) if not d["define"].endswith("_EXPORTS")]
        includes = [os.path.realpath(i["path"]) for i in group[0].get("includes", [])]

        # The new src/generated: the committed generated files CMake still
        # compiles at this ceiling, then whatever it generated into the build
        # tree (main set, Rys, MD), each keeping its path under its directory.
        old = os.path.realpath(os.path.join(tree, "src", "generated"))
        new = os.path.join(work, "generated")
        os.makedirs(new)
        placed = set()

        def put(src, rel):
            dst = os.path.join(new, rel)
            if os.path.exists(dst):
                sys.exit(f"make_fpm_tarball: {rel} would be written twice")
            os.makedirs(os.path.dirname(dst), exist_ok=True)
            shutil.copy2(src, dst)

        for p in glob.glob(os.path.join(old, "**", "*"), recursive=True):
            if not os.path.isfile(p):
                continue
            if os.path.realpath(p) in srcs:
                put(p, os.path.relpath(p, old))
                placed.add(os.path.realpath(p))
            elif not p.lower().endswith(".f90") and old in includes:
                put(p, os.path.relpath(p, old))  # the committed set's .inc files, when they are the ones included
        for gdir, sub in (("generated", ""), ("generated_rys", "rys"), ("generated_mdj", "mdj")):
            gpath = os.path.realpath(os.path.join(build, gdir))
            if not os.path.isdir(gpath):
                continue
            for p in glob.glob(os.path.join(gpath, "**", "*"), recursive=True):
                if not os.path.isfile(p):
                    continue
                rp = os.path.realpath(p)
                if rp in srcs:
                    put(p, os.path.join(sub, os.path.relpath(p, gpath)))
                    placed.add(rp)
                elif not p.lower().endswith(".f90") and gpath in includes:
                    put(p, os.path.join(sub, os.path.relpath(p, gpath)))
        missing = [p for p in srcs if not p.startswith(os.path.realpath(tree) + os.sep) and p not in placed]
        if missing:
            sys.exit("make_fpm_tarball: CMake compiles generated sources this did not place:\n  " +
                     "\n  ".join(sorted(missing)[:20]))
        shutil.rmtree(old)
        shutil.copytree(new, old)

        # Every source fpm will compile under src/ is one CMake compiles, or
        # behind TERCO_ENABLE_DLFIND (src/opt, as in the committed fpm build).
        for p in glob.glob(os.path.join(tree, "src", "**", "*.[Ff]90"), recursive=True):
            rel = os.path.relpath(p, tree)
            if rel.startswith(os.path.join("src", "generated")):
                continue
            if os.path.realpath(p) not in srcs and not rel.startswith(os.path.join("src", "opt")):
                sys.exit(f"make_fpm_tarball: fpm would compile {rel}, which CMake does not")

        # fpm.toml: CMake's definitions for this ceiling.
        toml = os.path.join(tree, "fpm.toml")
        text = open(toml).read()
        macros = "macros = [" + ", ".join(f'"{d}"' for d in defines) + "]"
        text, n = re.subn(r"^macros\s*=\s*\[[^\]]*\]", macros, text, flags=re.M)
        if n != 1:
            sys.exit("make_fpm_tarball: fpm.toml has no single macros line")
        text = text.replace("# fpm builds the SHIPPED configuration and only that one",
                            f"# THIS TARBALL is TERCO_LMAX={a.lmax}, generated already (FPM.md); in the\n"
                            "# repository, fpm builds the shipped configuration and only that one", 1)
        open(toml, "w").write(text)

        for d in DROP:
            p = os.path.join(tree, d)
            if os.path.isdir(p):
                shutil.rmtree(p)
            elif os.path.exists(p):
                os.remove(p)
        os.remove(os.path.join(tree, "tools", "fpm_generate.sh"))
        for f in DATA:
            shutil.copy2(os.path.join(root, "test", "basis_sets", f), os.path.join(tree, "basis_sets", f))

        nfiles = sum(len(fs) for _, _, fs in os.walk(os.path.join(tree, "src", "generated")))
        open(os.path.join(tree, "FPM.md"), "w").write(FPM_MD.format(
            version=version, lmax=a.lmax, ref=a.ref, defines=", ".join(defines), name=name, nfiles=nfiles,
            what={2: "s, p and d", 3: "s through f", 4: "s through g"}.get(a.lmax, f"l <= {a.lmax}")))

        path = os.path.join(out, name + ".tar.gz")
        with tarfile.open(path, "w:gz") as t:
            t.add(tree, arcname=name)
        print(f"make_fpm_tarball: {path} ({os.path.getsize(path) / 1e6:.1f} MB, {nfiles} generated files)")
    finally:
        if a.keep:
            print(f"make_fpm_tarball: work directory kept at {work}")
        else:
            shutil.rmtree(work, ignore_errors=True)


FPM_MD = """# terco {version}, TERCO_LMAX={lmax}, for fpm

This is the terco source at `{ref}` with every kernel already generated for
`TERCO_LMAX={lmax}` (basis functions {what}), the way a CMake configure at that
ceiling generates them: {nfiles} files in `src/generated/`. `fpm.toml`'s macros
are CMake's for it: {defines}. Nothing else needs generating; do not run the
repository's `tools/fpm_generate.sh` here (it is not included).

## Build and install

pic-mpi's manifest compiles MPI backends, so build with an MPI wrapper around
the Fortran compiler you want (`mpifort` below). The dependencies (pic,
pic-blas, pic-mpi, json-fortran, toml-f) are fetched by fpm from GitHub, so the
build needs network access.

GPU (NVHPC; `cc90` is an H100/H200, `cc80` an A100, `cc70` a V100), with the
SDK's own MPI wrapper (`$NVHPC_ROOT/comm_libs/mpi/bin/mpifort`) or a cluster
MPI built for nvfortran:

    mkdir -p blas && ln -sf $NVHPC_ROOT/compilers/lib/libblas.so blas/libopenblas.so
    fpm install --prefix $HOME/install/terco/{version} --compiler mpifort --profile release \\
      --flag "-stdpar=gpu -acc -gpu=cc90,mem:separate -DTRC_DEVICE" \\
      --link-flag "-L$PWD/blas -llapack -cudalib=cublas,cusolver"

The link flags are what CMake adds by itself. pic-blas's manifest links
`-lopenblas`, and the SDK ships OpenBLAS as `libblas.so` with the rest of
LAPACK in `liblapack.so`, hence the symlink and `-llapack` (skip the symlink
if an OpenBLAS is already on the linker path); `-cudalib` brings in cuBLAS and
cuSOLVER. At run time the MPI's own libraries must come first on
`LD_LIBRARY_PATH` (loading the SDK's or the cluster's MPI module does that).

CPU (gfortran 15 or newer; a correct and much slower library):

    fpm install --prefix $HOME/install/terco/{version} --compiler mpifort --profile release \\
      --flag "-ffree-line-length-none"

Measured on a V100 box (NVHPC 26.5, 4 parallel compiles): the LMAX=2
tarball takes 25 minutes; LMAX=3 and 4 take much longer, the per-class kernels
being many large files. The result matches the CMake build: gly50/cc-pVDZ
energies equal to 1e-11, Fock builds within 4%.

## Basis sets

`basis_sets/` holds the shipped library (Basis Set Exchange JSON, see its
README) plus the data files the minao and SAP guesses and Hirshfeld charges
read. fpm does not install data, so point the CLI at it:

    export TERCO_BASIS_PATH=$PWD/basis_sets
    terco water.xyz -b def2-svp

## Differences from the CMake build

CMake sets two compiler options on single files that fpm cannot express:

- ifx: `-fp-model=precise` on the translated XC functionals
  (`src/xc/functionals/`), without which a few of them lose accuracy. With
  ifx, add it to `--flag` for the whole build (correct, somewhat slower).
- nvfortran: `-Wptxas,--disable-optimizer-constants` on the same files, which
  keeps their literals out of the device constant bank. Not needed for the
  functionals in this release.

`--opt` (DL-FIND) is off, as in a default CMake build; README.md, "fpm",
explains how to turn it on.
"""

if __name__ == "__main__":
    main()
