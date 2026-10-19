#!/usr/bin/env python3
#
# The `terco` command over MPI ranks: run on one rank and two, and check
# that the two agree -- the same thing scf_mpi.F90 checks for the library's
# own Fortran and C entries, here for the thing a user actually runs.
#
# Water RHF and B3LYP must match between np 1 and np 2 to 1e-8 (scf_mpi's
# own tolerance); wB97X-V (VV10) must be refused cleanly, on every rank, at
# np 2. RI-MP2's correlation energy must match between np 1 and np 2 to
# 1e-10, with B distributed over the ranks and passed round them, both with
# the build's default transport and staged through the host
# (TERCO_MPI_DEVICE=0).
#
import re
import subprocess
import sys

TOL = 1.0e-8
ENERGY_RE = re.compile(r"SCF ENERGY\s+(-?\d+\.\d+)")
MP2_RE = re.compile(r"MP2 CORRELATION\s+(-?\d+\.\d+)")
MP2_TOL = 1.0e-10


def run(mpiexec, npflag, np, terco_cli, xyz, basis_dir, *args, extra_env=None):
    import os
    cmd = [mpiexec, npflag, str(np), terco_cli, xyz, "-b", "6-31g", *args]
    env = {**os.environ, "TERCO_BASIS_PATH": basis_dir, **(extra_env or {})}
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
            print(f"cli_mpi_check: {label} -- could not find SCF ENERGY in output")
            ok = False
            continue
        # np 2 must print the result once, not twice.
        if p2.stdout.count("SCF ENERGY") != 1:
            print(f"cli_mpi_check: {label} -- np 2 printed SCF ENERGY "
                  f"{p2.stdout.count('SCF ENERGY')} times, not once")
            ok = False
        diff = abs(e1 - e2)
        print(f"cli_mpi_check: {label}  np1 = {e1:.12f}  np2 = {e2:.12f}  diff = {diff:.2e}")
        if diff > TOL:
            ok = False

    # RI-MP2 with B split over the ranks: the default transport, then
    # through the host. OpenMPI passes the environment to every rank.
    m1 = run(mpiexec, npflag, 1, terco_cli, xyz, basis_dir, "--rimp2")
    mm1 = MP2_RE.search(m1.stdout)
    for label, env in [("default", None), ("host-staged", {"TERCO_MPI_DEVICE": "0"})]:
        m2 = run(mpiexec, npflag, 2, terco_cli, xyz, basis_dir, "--rimp2", extra_env=env)
        mm2 = MP2_RE.search(m2.stdout)
        if m1.returncode != 0 or m2.returncode != 0 or mm1 is None or mm2 is None:
            print(f"cli_mpi_check: RI-MP2 ({label}) did not run cleanly "
                  f"(np1 rc={m1.returncode}, np2 rc={m2.returncode})")
            print(m1.stdout, m1.stderr, m2.stdout, m2.stderr)
            ok = False
            continue
        c1, c2 = float(mm1.group(1)), float(mm2.group(1))
        diff = abs(c1 - c2)
        print(f"cli_mpi_check: RI-MP2 ({label})  np1 = {c1:.12f}  np2 = {c2:.12f}  diff = {diff:.2e}")
        if diff > MP2_TOL:
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
