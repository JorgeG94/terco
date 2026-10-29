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
# MULTI-GPU header), not just the evenly-split np 2 case -- measured,
# neither grid's point count is a multiple of 3 (see the np3 block below
# for the actual counts), so this is not an assumption.
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
# Hirshfeld-condensed Fukui indices (three SCFs on one handle, each rank
# running them collectively, the Hirshfeld analysis replicated) must match
# between np 1, np 2 and np 3 to 1e-8, read from --json, and the table must
# be printed once.
#
# RI-MP2's correlation energy must match between np 1 and np 2 to 1e-10,
# with (mu nu|P) split by auxiliary shell and B distributed over the ranks
# and passed round them, both with the build's default transport and
# staged through the host (TERCO_MPI_DEVICE=0), and between np 1 and np 3
# with uneven occupied blocks and several auxiliary blocks per rank.
#
import json
import os
import re
import subprocess
import sys
import tempfile

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
    mpiexec, npflag, terco_cli, xyz, basis_dir, nh2_xyz = sys.argv[1:7]
    ok = True

    # check_exc/check_np3 are explicit per case rather than keyed off a
    # label substring: both are VV10 properties (the exchange-correlation
    # row only has a meaningful E_nlc contribution, and only VV10's own
    # grid is split, when wB97X-V is running), not something RHF/B3LYP
    # need or NH2 happens to share a name with. np3 is deliberately run
    # only for the two VV10 cases: RHF/B3LYP have no grid to split
    # unevenly, so a remainder check there would exercise nothing new.
    cases = [
        ("RHF", xyz, [], False, False),
        ("B3LYP", xyz, ["-f", "b3lyp"], False, False),
        ("wB97X-V", xyz, ["-f", "wb97x-v"], True, True),
        ("wB97X-V (UKS, NH2)", nh2_xyz, ["-f", "wb97x-v"], True, True),
    ]
    for label, mol, extra, check_exc, check_np3 in cases:
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

        if check_exc:
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

        if check_np3:
            # np 3: the grid's point count is not a multiple of 3 (water's
            # own VV10 grid: 10124 points, 10124 % 3 == 2; NH2's: 10064
            # points, 10064 % 3 == 2 -- both measured at this program's
            # default --nlc-grid level 1), so the "first (points mod
            # ranks) ranks get one extra point" remainder handling runs
            # too, not just the evenly-split np 2 case above.
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

    # Fukui over ranks: wB97X-V (VV10's rank split runs in all three SCFs)
    # from the minao guess, which the ions inherit.
    fukui = {}
    with tempfile.TemporaryDirectory() as tmp:
        for np in (1, 2, 3):
            js = os.path.join(tmp, f"fukui{np}.json")
            pf = run(mpiexec, npflag, np, terco_cli, xyz, basis_dir, "-f", "wb97x-v",
                     "--guess", "minao", "--fukui", "hirshfeld", "--json", js)
            if pf.returncode != 0 or not os.path.exists(js):
                print(f"cli_mpi_check: Fukui np{np} did not run cleanly (rc={pf.returncode})")
                print(pf.stdout, pf.stderr)
                ok = False
                continue
            if pf.stdout.count("f+ is where the molecule accepts charge") != 1:
                print(f"cli_mpi_check: Fukui np{np} printed its table "
                      f"{pf.stdout.count('f+ is where the molecule accepts charge')} times, not once")
                ok = False
            with open(js) as fh:
                fukui[np] = json.load(fh)["fukui"]
    if 1 in fukui:
        for np in (2, 3):
            if np not in fukui:
                continue
            worst = max(abs(a - b) for key in ("f+", "f-")
                        for a, b in zip(fukui[1][key], fukui[np][key]))
            print(f"cli_mpi_check: Fukui (wB97X-V, minao)  np1 vs np{np}  worst |df| = {worst:.2e}")
            if worst > TOL:
                ok = False

    # RI-MP2 with (mu nu|P) split by auxiliary shell and B by occupied
    # block: the default transport, then through the host, then np 3 with
    # small auxiliary blocks, so the occupied blocks are uneven (water
    # 6-31G correlates 4 orbitals: 2, 1, 1) and each rank's auxiliary
    # shells take several blocks. OpenMPI passes the environment to every
    # rank.
    m1 = run(mpiexec, npflag, 1, terco_cli, xyz, basis_dir, "--rimp2")
    mm1 = MP2_RE.search(m1.stdout)
    for label, np, env, extra in [("default", 2, None, []),
                                  ("host-staged", 2, {"TERCO_MPI_DEVICE": "0"}, []),
                                  ("np 3, aux blocks of 8", 3, None, ["--aux-block", "8"])]:
        m2 = run(mpiexec, npflag, np, terco_cli, xyz, basis_dir, "--rimp2", *extra, extra_env=env)
        mm2 = MP2_RE.search(m2.stdout)
        if m1.returncode != 0 or m2.returncode != 0 or mm1 is None or mm2 is None:
            print(f"cli_mpi_check: RI-MP2 ({label}) did not run cleanly "
                  f"(np1 rc={m1.returncode}, np{np} rc={m2.returncode})")
            print(m1.stdout, m1.stderr, m2.stdout, m2.stderr)
            ok = False
            continue
        c1, c2 = float(mm1.group(1)), float(mm2.group(1))
        diff = abs(c1 - c2)
        print(f"cli_mpi_check: RI-MP2 ({label})  np1 = {c1:.12f}  np{np} = {c2:.12f}  diff = {diff:.2e}")
        if diff > MP2_TOL:
            ok = False

    print("cli_mpi_check: PASS" if ok else "cli_mpi_check: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
