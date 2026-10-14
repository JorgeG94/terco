#!/usr/bin/env python3
#
# The `terco` command over MPI ranks: run on one rank and two, and check
# that the two agree -- the same thing scf_mpi.F90 checks for the library's
# own Fortran and C entries, here for the thing a user actually runs.
#
# Water RHF and B3LYP must match between np 1 and np 2 to 1e-8 (scf_mpi's
# own tolerance); wB97X-V (VV10) must be refused cleanly, on every rank, at
# np 2.
#
import re
import subprocess
import sys

TOL = 1.0e-8
ENERGY_RE = re.compile(r"TOTAL ENERGY\s+(-?\d+\.\d+)")


def run(mpiexec, npflag, np, terco_cli, xyz, basis_dir, *args):
    import os
    cmd = [mpiexec, npflag, str(np), terco_cli, xyz, "-b", "6-31g", *args]
    env = {**os.environ, "TERCO_BASIS_PATH": basis_dir}
    return subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300)


def energy(out):
    m = ENERGY_RE.search(out)
    return float(m.group(1)) if m else None


def main():
    mpiexec, npflag, terco_cli, xyz, basis_dir = sys.argv[1:6]
    ok = True

    for label, extra in [("RHF", []), ("B3LYP", ["-f", "b3lyp"])]:
        p1 = run(mpiexec, npflag, 1, terco_cli, xyz, basis_dir, *extra)
        p2 = run(mpiexec, npflag, 2, terco_cli, xyz, basis_dir, *extra)
        if p1.returncode != 0 or p2.returncode != 0:
            print(f"cli_mpi_check: {label} did not run cleanly "
                  f"(np1 rc={p1.returncode}, np2 rc={p2.returncode})")
            print(p1.stdout, p1.stderr, p2.stdout, p2.stderr)
            ok = False
            continue
        e1, e2 = energy(p1.stdout), energy(p2.stdout)
        if e1 is None or e2 is None:
            print(f"cli_mpi_check: {label} -- could not find TOTAL ENERGY in output")
            ok = False
            continue
        # np 2 must print the result once, not twice.
        if p2.stdout.count("TOTAL ENERGY") != 1:
            print(f"cli_mpi_check: {label} -- np 2 printed TOTAL ENERGY "
                  f"{p2.stdout.count('TOTAL ENERGY')} times, not once")
            ok = False
        diff = abs(e1 - e2)
        print(f"cli_mpi_check: {label}  np1 = {e1:.12f}  np2 = {e2:.12f}  diff = {diff:.2e}")
        if diff > TOL:
            ok = False

    # wB97X-V needs VV10, which the library refuses over more than one rank
    # -- every rank must say so and exit non-zero, not hang.
    pv = run(mpiexec, npflag, 2, terco_cli, xyz, basis_dir, "-f", "wb97x-v")
    refused = pv.returncode != 0 and "VV10" in (pv.stdout + pv.stderr)
    print(f"cli_mpi_check: wB97X-V at np 2  rc={pv.returncode}  refused={refused}")
    if not refused:
        print(pv.stdout, pv.stderr)
        ok = False

    print("cli_mpi_check: PASS" if ok else "cli_mpi_check: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
