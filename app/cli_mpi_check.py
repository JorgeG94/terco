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
# between np 1, np 2 and np 3 to 1e-8, read from --json, and the table and
# each SCF stage label must be printed once.
#
# RI-MP2's correlation energy must match between np 1 and np 2 to 1e-10,
# with (mu nu|P) split by auxiliary shell and B distributed over the ranks
# and passed round them, both with the build's default transport and
# staged through the host (TERCO_MPI_DEVICE=0), and between np 1 and np 3
# with uneven occupied blocks and several auxiliary blocks per rank.
# RI-CCSD(T) (water, cc-pVDZ with cc-pVDZ-RIFIT, frozen core 1): the CCSD
# correlation energy and E(T) read from --json must match between np 1, np 2
# and np 3 to CC_TOL (measured 7e-15; the check prints the spread it finds), with the ladder of every iteration and the (T) shared. What
# only a CLI run over ranks shows is checked too: rank 0's energy box says
# "ladder and (T) shared by N ranks" (and a one-rank run does not), the step
# table carries "rank 0 of N", the result is printed once; and the options
# refused over ranks (an explicit reference or tiled (T) routine, rejected by
# the CLI before any work; an explicit assembled ladder, refused by the
# library after the SCF on every rank) end the whole run with a non-zero exit
# status and a reason, within the timeout, instead of leaving a rank waiting.
#
import json
import os
import re
import signal
import subprocess
import sys
import tempfile

TOL = 1.0e-8
ENERGY_RE = re.compile(r"TOTAL ENERGY\s+(-?\d+\.\d+)")
EXC_RE = re.compile(r"exchange-correlation\s+(-?\d+\.\d+)")
MP2_RE = re.compile(r"RI-MP2 CORRELATION\s+(-?\d+\.\d+)")
MP2_TOL = 1.0e-10
CC_TOL = 1.0e-12
# RI-CCSD correlation energy and E(T) between rank counts. Measured with these settings (nvfortran 26.5 host build,
# water cc-pVDZ / cc-pVDZ-RIFIT, np 2 and 3 against np 1): 7.0e-15 Eh for E_CCSD and 3.8e-16 Eh for E(T); the check
# prints what it finds (the "diff E_CCSD" and "diff E(T)" columns). The settings below converge the amplitudes to
# 1e-10, the energy to 1e-12 and the SCF to 1e-13 / 1e-11; what is left is the multi-rank SCF (about 1e-14 Eh,
# docs/sphinx/cli.md) and the order of the ladder's sum. 1e-12 leaves two orders for another MPI or BLAS.
CC_ARGS = ["--frozen", "1", "--econv", "1e-13", "--gconv", "1e-11", "--cc-econv", "1e-12", "--cc-tconv", "1e-10"]


def run(mpiexec, npflag, np, terco_cli, xyz, basis_dir, *args, extra_env=None):
    cmd = [mpiexec, npflag, str(np), terco_cli, xyz, "-b", "6-31g", *args]
    env = {**os.environ, "TERCO_BASIS_PATH": basis_dir, **(extra_env or {})}
    return subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=300)


def _text(x):
    """subprocess puts bytes (or None) in TimeoutExpired.stdout/.stderr even with text=True"""
    if x is None:
        return ""
    return x.decode(errors="replace") if isinstance(x, bytes) else x


def run_group(cmd, env, timeout):
    """Run `cmd` in its own session; on a timeout kill the whole group (no orphan ranks) and return rc -999 with
    whatever the ranks had printed."""
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env,
                            start_new_session=True)
    try:
        out, err = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired as exc:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        out, err = proc.communicate()
        out = _text(exc.stdout) + _text(out)
        err = _text(exc.stderr) + _text(err) + "\nTIMEOUT (a rank is waiting?)"
        return subprocess.CompletedProcess(cmd, -999, out, err)
    return subprocess.CompletedProcess(cmd, proc.returncode, out, err)


def run_cc(mpiexec, npflag, np, terco_cli, xyz, basis_dir, *args, timeout=300):
    cmd = [mpiexec, npflag, str(np), terco_cli, xyz, "-b", os.path.join(basis_dir, "cc-pvdz.json"),
           "--aux", os.path.join(basis_dir, "cc-pvdz-rifit.json"), *CC_ARGS, *args]
    env = {**os.environ, "TERCO_BASIS_PATH": basis_dir}
    return run_group(cmd, env, timeout)


def energy(out):
    m = ENERGY_RE.search(out)
    return float(m.group(1)) if m else None


def exc_energy(out):
    m = EXC_RE.search(out)
    return float(m.group(1)) if m else None


def check_riccsd(mpiexec, npflag, terco_cli, xyz, basis_dir):
    ok = True
    res = {}
    with tempfile.TemporaryDirectory() as tmp:
        for np in (1, 2, 3):
            js = os.path.join(tmp, f"cc{np}.json")
            p = run_cc(mpiexec, npflag, np, terco_cli, xyz, basis_dir, "--riccsd-t", "--json", js)
            if p.returncode != 0 or not os.path.exists(js):
                print(f"cli_mpi_check: RI-CCSD(T) np{np} did not run cleanly (rc={p.returncode})")
                print(p.stdout, p.stderr)
                ok = False
                continue
            with open(js) as fh:
                d = json.load(fh)
            res[np] = (d["ri-ccsd correlation energy"], d["ri-ccsd triples correction"])
            if p.stdout.count("TOTAL RI-CCSD(T) ENERGY") != 1:
                print(f"cli_mpi_check: RI-CCSD(T) np{np} printed the total "
                      f"{p.stdout.count('TOTAL RI-CCSD(T) ENERGY')} times, not once")
                ok = False
            shared = f"ladder and (T) shared by {np} ranks"
            has_note = shared in p.stdout
            if (np > 1) != has_note:
                print(f"cli_mpi_check: RI-CCSD(T) np{np}: the note '{shared}' is "
                      f"{'missing' if np > 1 else 'printed on one rank'}")
                ok = False
            table = f"rank 0 of {np}"
            if table not in p.stdout:
                print(f"cli_mpi_check: RI-CCSD(T) np{np}: the step table does not say '{table}'")
                ok = False
            if np > 1 and "amplitudes to ranks" not in p.stdout:
                print(f"cli_mpi_check: RI-CCSD(T) np{np}: the step table has no transfer rows")
                ok = False
    if 1 in res:
        for np in (2, 3):
            if np not in res:
                continue
            dc, dt = abs(res[1][0] - res[np][0]), abs(res[1][1] - res[np][1])
            print(f"cli_mpi_check: RI-CCSD(T)  np1 E_CCSD = {res[1][0]:.12f}  E(T) = {res[1][1]:.12f}  "
                  f"np{np} diff E_CCSD = {dc:.2e}  diff E(T) = {dt:.2e}")
            if dc > CC_TOL or dt > CC_TOL:
                ok = False

    # Refusals over ranks: non-zero exit, a reason, no hang.
    # each matched against a distinctive fragment of the refusal itself (app/trc_run_cc.F90 for the first two,
    # trc_rccsd_driver's ladder_for_ranks for the third), not against a word normal output also contains
    refusals = [
        ("explicit reference (T) at np 2", 2, ["--riccsd-t", "--cc-triples", "reference"],
         "reference and tiled run on one rank"),
        ("explicit tiled (T) at np 3", 3, ["--riccsd-t", "--cc-triples", "tiled"],
         "reference and tiled run on one rank"),
        ("explicit assembled ladder at np 2", 2, ["--riccsd-t", "--cc-ladder", "assembled"],
         "needs the half form"),
    ]
    for label, np, args, word in refusals:
        p = run_cc(mpiexec, npflag, np, terco_cli, xyz, basis_dir, *args, timeout=120)
        text = p.stdout + p.stderr
        good = p.returncode not in (0, -999) and word in text and "TOTAL RI-CCSD" not in p.stdout
        print(f"cli_mpi_check: refusal ({label})  rc = {p.returncode}  "
              f"{'ok' if good else 'NOT REFUSED AS EXPECTED'}")
        if not good:
            print(p.stdout, p.stderr)
            ok = False
    return ok


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
            # Each stage label once too: they are logged outside the SCF,
            # so nothing else gates them to rank 0.
            for stage in ("N (reference", "N+1 (anion)", "N-1 (cation)"):
                count = pf.stdout.count("Fukui: " + stage)
                if count != 1:
                    print(f"cli_mpi_check: Fukui np{np} printed the '{stage}' stage label "
                          f"{count} times, not once")
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

    ok = check_riccsd(mpiexec, npflag, terco_cli, xyz, basis_dir) and ok

    print("cli_mpi_check: PASS" if ok else "cli_mpi_check: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
