#!/usr/bin/env python3
#
# The `terco` command over MPI ranks: run on one rank and more, and check
# that they agree -- the same thing scf_mpi.F90 checks for the library's
# own Fortran and C entries, here for the thing a user actually runs.
#
# Water RHF, B3LYP and wB97X-V (RKS, VV10 split by a contiguous range of
# pass 2's grid points per rank, one allreduce) and NH2 radical wB97X-V
# (UKS, same split) must all match between np 1, np 2 and np 3 to 1e-8
# (scf_mpi's own tolerance) -- both the total energy and, for wB97X-V, the
# exchange-correlation energy the "exchange-correlation" row prints
# (semilocal XC plus VV10's E_nlc, the quantity pass 2's rank split and
# allreduce actually touch). np 3 exercises the "first (points mod ranks)
# ranks get one extra point" remainder handling (see trc_vv10.F90's own
# MULTI-GPU header), not just the evenly-split np 2 case -- both grids'
# point counts are not multiples of 3.
#
# NH2, not OH, is the open-shell case: OH's SOMO is a doubly-degenerate pi
# orbital, so each SCF run is free to converge to an arbitrary mixture of
# the two pi occupations, and a run-to-run perturbation as small as grid-
# batch ordering is enough to pick a different one -- measured, repeated
# np 1 runs of OH wB97X-V alone (same rank count, same binary, same
# machine) disagreed by up to ~2e-7 Eh, and tightening --gconv/--econv by
# two orders of magnitude did not shrink it, so it was never about
# convergence or about MPI rank count. NH2's SOMO is the non-degenerate
# b1 orbital (see test/nh2.xyz's own comment), and is bit-identical run-
# to-run and across np 1/2/3 at the same 1e-8 TOL every RKS case uses.
#
# RI-MP2's correlation energy must match between np 1 and np 2 to 1e-10,
# with B distributed over the ranks and passed round them, both with the
# build's default transport and staged through the host
# (TERCO_MPI_DEVICE=0).
#
import os
import re
import subprocess
import sys

TOL = 1.0e-8
ENERGY_RE = re.compile(r"TOTAL ENERGY\s+(-?\d+\.\d+)")
EXC_RE = re.compile(r"exchange-correlation\s+(-?\d+\.\d+)")
MP2_RE = re.compile(r"MP2 CORRELATION\s+(-?\d+\.\d+)")
MP2_TOL = 1.0e-10


def run(mpiexec, npflag, np, terco_cli, xyz, basis_dir, *args, extra_env=None):
    cmd = [mpiexec, npflag, str(np), terco_cli, xyz, "-b", "6-31g", *args]
    env = {**os.environ, "TERCO_BASIS_PATH": basis_dir, **(extra_env or {})}
    return subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300)


def energy(out):
    m = ENERGY_RE.search(out)
    return float(m.group(1)) if m else None


def exc_energy(out):
    m = EXC_RE.search(out)
    return float(m.group(1)) if m else None


def main():
    mpiexec, npflag, terco_cli, xyz, basis_dir = sys.argv[1:6]
    nh2_xyz = os.path.join(os.path.dirname(xyz), "nh2.xyz")
    ok = True

    cases = [
        ("RHF", xyz, []),
        ("B3LYP", xyz, ["-f", "b3lyp"]),
        ("wB97X-V", xyz, ["-f", "wb97x-v"]),
        ("wB97X-V (UKS, NH2)", nh2_xyz, ["-f", "wb97x-v"]),
    ]
    for label, mol, extra in cases:
        p1 = run(mpiexec, npflag, 1, terco_cli, mol, basis_dir, *extra)
        p2 = run(mpiexec, npflag, 2, terco_cli, mol, basis_dir, *extra)
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

        if "wB97X-V" in label:
            # The exchange-correlation row carries VV10's E_nlc -- the
            # quantity pass 2's rank split and allreduce actually touch --
            # separately from the total, which also depends on H and G.
            x1, x2 = exc_energy(p1.stdout), exc_energy(p2.stdout)
            if x1 is None or x2 is None:
                print(f"cli_mpi_check: {label} -- could not find exchange-correlation in output")
                ok = False
                continue
            xdiff = abs(x1 - x2)
            print(f"cli_mpi_check: {label} exchange-correlation (incl. VV10's E_nlc)  "
                  f"np1 = {x1:.12f}  np2 = {x2:.12f}  diff = {xdiff:.2e}")
            if xdiff > TOL:
                ok = False

            # np 3: the grid's point count is not a multiple of 3, so the
            # "first (points mod ranks) ranks get one extra point" remainder
            # handling runs too, not just the evenly-split np 2 case above.
            p3 = run(mpiexec, npflag, 3, terco_cli, mol, basis_dir, *extra)
            if p3.returncode != 0:
                print(f"cli_mpi_check: {label} np3 did not run cleanly (rc={p3.returncode})")
                print(p3.stdout, p3.stderr)
                ok = False
                continue
            e3, x3 = energy(p3.stdout), exc_energy(p3.stdout)
            if e3 is None or x3 is None:
                print(f"cli_mpi_check: {label} np3 -- could not find TOTAL ENERGY/"
                      f"exchange-correlation in output")
                ok = False
                continue
            diff3, xdiff3 = abs(e1 - e3), abs(x1 - x3)
            print(f"cli_mpi_check: {label} np3  total diff = {diff3:.2e}  "
                  f"exchange-correlation diff = {xdiff3:.2e}")
            if diff3 > TOL or xdiff3 > TOL:
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

    print("cli_mpi_check: PASS" if ok else "cli_mpi_check: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
