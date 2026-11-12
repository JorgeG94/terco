#!/usr/bin/env python3
"""Performance-regression bench for terco: terco against ITSELF over time.

Runs the cases of cases.toml through a built `terco` binary with --json,
takes the median of every stage over --repeats runs (after --warmup discarded
runs) and compares with the golden CSV of the current GPU architecture
(goldens/<arch>.csv). Only the stages a case lists under `judge` can FAIL: a
stage slower than its golden by more than its fail ratio. Every other recorded
stage is reported as INFO. An energy that moved by more than `energy_tol` is a
FAIL whatever the timings say, so a faster wrong answer fails.

Exit status: 0 no regression, 1 at least one FAIL, 2 usage error or at least
one case that errored, did not converge or was refused (2 wins over 1: the
bench did not finish its job).

Standard library only (Python >= 3.11 for tomllib). See README.md.
"""

import argparse
import csv
import datetime
import json
import os
import platform
import re
import signal
import socket
import statistics
import subprocess
import sys
import tempfile

try:
    import tomllib
except ImportError:  # pragma: no cover
    tomllib = None

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))

PROVENANCE = ["binary_commit", "checkout_commit", "date", "host", "cpu_model", "gpu_name", "driver", "cuda",
              "omp_threads", "mkl_threads", "arch", "build", "compiler"]
GOLDEN_FIELDS = ["case", "stage", "median_s", "min_s", "max_s", "repeats", "energy_hartree", "iterations",
                 "gradient_rms"] + PROVENANCE

DEFAULTS = {
    "fail": 1.10,
    "warn": 1.05,
    "min_seconds": 0.05,
    "energy_tol": 1.0e-8,
    "gradient_tol": 1.0e-6,
    "max_seconds": 300,
}

SPREAD_REPEATS = 5     # a golden recorded with at least this many repeats carries a usable spread
SPREAD_FACTOR = 4.0    # fail at 1 + SPREAD_FACTOR * (max - min) / median, if that is above the configured ratio

# Stages computed from others, per repeat: name -> (numerator, denominator).
DERIVED = {"scf.fock_per_iteration": ("scf.fock_total", "scf.n_fock_builds")}

# The `terco --version` build facts a golden must agree with (gpu_arch only for device builds).
BUILD_KEYS = ("device", "gpu_arch", "lmax", "df_lmax", "df_lmax_aux", "cart_order", "mpi", "multicore")

ERROR_VERDICTS = ("ERROR", "NOT CONVERGED", "BUILD MISMATCH")


class UsageError(Exception):
    """A problem with the invocation or the environment: exit status 2."""


class RunError(Exception):
    """A terco run that failed, timed out, did not converge or wrote no JSON."""

    def __init__(self, message, verdict="ERROR", log=None):
        super().__init__(message)
        self.verdict = verdict
        self.log = log


# ----------------------------------------------------------------------------
# configuration
# ----------------------------------------------------------------------------

def load_config(path):
    """cases.toml -> (defaults, per-stage thresholds, list of case dicts)."""
    if tomllib is None:
        raise UsageError("Python >= 3.11 is needed for tomllib")
    try:
        with open(path, "rb") as fh:
            raw = tomllib.load(fh)
    except (OSError, tomllib.TOMLDecodeError) as exc:
        raise UsageError(f"cannot read {path}: {exc}")
    defaults = dict(DEFAULTS)
    defaults.update(raw.get("defaults", {}))
    thresholds = raw.get("thresholds", {})
    cases = raw.get("case", [])
    base = os.path.dirname(os.path.abspath(path))
    seen = set()
    for case in cases:
        for key in ("name", "geometry", "args", "stages"):
            if key not in case:
                raise UsageError(f"{path}: a case lacks '{key}': {case}")
        if case["name"] in seen:
            raise UsageError(f"{path}: duplicate case name {case['name']}")
        seen.add(case["name"])
        case["geometry_path"] = os.path.join(base, case["geometry"])
        case.setdefault("needs", [])
        case.setdefault("quick", False)
        case.setdefault("judge", [])
        case.setdefault("iteration_verdict", "WARN")
        stray = [s for s in case["judge"] if s not in case["stages"]]
        if stray:
            raise UsageError(f"{path}: case {case['name']} judges {stray}, which are not in its stages")
        if case["iteration_verdict"] not in ("WARN", "FAIL"):
            raise UsageError(f"{path}: case {case['name']}: iteration_verdict must be WARN or FAIL")
    return defaults, thresholds, cases


def thresholds_for(defaults, per_stage, case, stage):
    """fail/warn/min_seconds for one (case, stage): case override, then the
    per-stage table, then the global default."""
    out = {k: defaults[k] for k in ("fail", "warn", "min_seconds")}
    # A case may also set fail, warn or min_seconds for all of its stages (the
    # water-hexamer DF cases lower the floor: their stages are milliseconds).
    for src in (per_stage.get(stage, {}), case, case.get("thresholds", {}).get(stage, {})):
        for k in out:
            if k in src:
                out[k] = src[k]
    return out


def effective_ratios(thr, grow):
    """(fail, warn, source) for one stage. A golden row recorded with at least
    SPREAD_REPEATS repeats gives fail = max(configured fail, 1 + 4 x (max - min) / median)
    and warn half-way between 1 and that; otherwise the configured ratios apply."""
    fail, warn = thr["fail"], thr["warn"]
    if grow is not None and (grow.get("repeats") or 0) >= SPREAD_REPEATS and (grow.get("median_s") or 0) > 0:
        spread = (grow["max_s"] - grow["min_s"]) / grow["median_s"]
        fail = max(fail, 1.0 + SPREAD_FACTOR * spread)
        return fail, 1.0 + (fail - 1.0) / 2.0, "spread"
    return fail, warn, "configured"


# ----------------------------------------------------------------------------
# the binary: build facts, GPU, architecture slug
# ----------------------------------------------------------------------------

def discover_terco(explicit):
    """The terco binary: --terco, or a build*/app/terco under the repo."""
    if explicit:
        if not (os.path.isfile(explicit) and os.access(explicit, os.X_OK)):
            raise UsageError(f"--terco {explicit}: not an executable file")
        return os.path.abspath(explicit)
    found = []
    try:
        names = sorted(os.listdir(REPO))
    except OSError:
        names = []
    for name in names:
        cand = os.path.join(REPO, name, "app", "terco")
        if name.startswith("build") and os.path.isfile(cand) and os.access(cand, os.X_OK):
            found.append(cand)
    if not found:
        raise UsageError(
            f"no terco binary found: looked for {REPO}/build*/app/terco. Build terco (README.md) or pass --terco PATH")
    if len(found) > 1:
        newest = max(found, key=os.path.getmtime)
        print(f"note: several builds found ({', '.join(found)}); using the newest, {newest}. "
              f"Pass --terco to choose.", file=sys.stderr)
        return newest
    return found[0]


def parse_build_line(text):
    """The 'build: key=value ...' line of `terco --version` -> dict, or None.
    compiler is two words ('GNU 15.1.0'), so it takes the rest of the line."""
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("build:"):
            body = line[len("build:"):].strip()
            info = {}
            m = re.search(r"\bcompiler=(.*)$", body)
            if m:
                info["compiler"] = m.group(1).strip()
                body = body[:m.start()]
            for tok in body.split():
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    info[k] = v
            return info
    return None


def probe_build(terco):
    try:
        proc = subprocess.run([terco, "--version"], capture_output=True, text=True, timeout=60)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise UsageError(f"cannot run {terco} --version: {exc}")
    info = parse_build_line(proc.stdout + proc.stderr)
    if info is None:
        raise UsageError(f"{terco} --version printed no 'build:' line: this terco is too old for the bench "
                         f"(it has no \"timings\" in --json either). Rebuild from a tree with them.")
    info["device_word"] = info.get("device", "no")
    info["device"] = info["device_word"] == "yes"
    try:
        info["lmax_int"] = int(info.get("lmax", "0"))
    except ValueError:
        info["lmax_int"] = 0
    return info


def build_capabilities(info):
    caps = set()
    if info["device"]:
        caps.add("device")
    if info["lmax_int"] >= 3:
        caps.add("lmax3")
    return caps


def build_string(info):
    """The kernel-shaping facts of a build as 'k=v;k=v', for the golden's build column."""
    parts = []
    for key in BUILD_KEYS:
        if key == "device":
            parts.append(f"device={info['device_word']}")
        elif key == "gpu_arch" and not info["device"]:
            continue
        elif key in info:
            parts.append(f"{key}={info[key]}")
    return ";".join(parts)


def parse_build_string(text):
    out = {}
    for tok in (text or "").split(";"):
        if "=" in tok:
            k, v = tok.split("=", 1)
            out[k] = v
    return out


def compiler_family(text):
    return (text or "").split()[0].lower() if (text or "").strip() else ""


def build_mismatch(golden_build, golden_compiler, now_build, now_compiler):
    """What differs between a golden's build and the current binary's, as
    sentences. Keys either side does not carry (an older golden) are not compared."""
    gold, now = parse_build_string(golden_build), parse_build_string(now_build)
    out = []
    for key in BUILD_KEYS:
        if key in gold and key in now and gold[key] != now[key]:
            out.append(f"{key}: golden {gold[key]}, this binary {now[key]}")
    gf, nf = compiler_family(golden_compiler), compiler_family(now_compiler)
    if gf and nf and gf != nf:
        out.append(f"compiler family: golden {golden_compiler}, this binary {now_compiler}")
    return out


def resolve_gpu(device, gpu_arg, allow_gpu0, environ):
    """The CUDA_VISIBLE_DEVICES value the child gets, or None for a host build.
    Never GPU 0 unless --allow-gpu0: this box's ctest pins GPU 0 and other
    work runs there. A GPU-<uuid> (or MIG-) value from the environment is
    allowed, but it cannot be checked against GPU 0 (see gpu_is_uuid).
    Raises UsageError with the reason."""
    if not device:
        return None
    if gpu_arg is not None:
        if gpu_arg < 0:
            raise UsageError(f"--gpu {gpu_arg}: a GPU index cannot be negative")
        value = str(gpu_arg)
    elif environ.get("CUDA_VISIBLE_DEVICES", "").strip():
        value = environ["CUDA_VISIBLE_DEVICES"].strip()
    else:
        raise UsageError("this is a GPU build and no GPU was chosen: pass --gpu N (never 0 unless you also pass "
                         "--allow-gpu0) or set CUDA_VISIBLE_DEVICES. The bench does not pick a GPU on its own.")
    first = value.split(",")[0].strip()
    if first.startswith(("GPU-", "MIG-")):
        return value
    if not first.isdigit():
        raise UsageError(f"CUDA_VISIBLE_DEVICES={value!r} is neither GPU indices nor GPU-<uuid>: refusing to guess")
    if first == "0" and not allow_gpu0:
        raise UsageError("refusing to run on GPU 0 (shared with ctest and other work): choose another GPU, "
                         "or pass --allow-gpu0 to run there deliberately")
    return value


def gpu_is_uuid(value):
    return value is not None and value.split(",")[0].strip().startswith(("GPU-", "MIG-"))


def slugify(text):
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")


def gpu_slug(name):
    """Architecture slug (the family: v100, a100, h100, h200, gh200) from the name nvidia-smi reports."""
    low = name.lower()
    if re.search(r"gh200|grace[ -]?hopper", low):
        return "gh200"
    for tag in ("h200", "h100", "a100", "v100"):
        if re.search(tag, low):
            return tag
    cleaned = re.sub(r"\bnvidia\b|\btesla\b|\bgeforce\b", " ", low)
    return slugify(cleaned) or "unknown-gpu"


def cpu_model_name(cpuinfo_text):
    for line in cpuinfo_text.splitlines():
        if line.lower().startswith("model name"):
            return line.split(":", 1)[1].strip()
    return ""


def cpu_slug(cpuinfo_text):
    """cpu-<model slug> from the 'model name' line of /proc/cpuinfo."""
    model = cpu_model_name(cpuinfo_text)
    model = re.sub(r"\((r|tm)\)", " ", model, flags=re.I)
    model = re.sub(r"@.*$", " ", model)
    model = re.sub(r"\b(cpu|processor)\b", " ", model, flags=re.I)
    return "cpu-" + (slugify(model) or "unknown")


def read_cpuinfo():
    try:
        with open("/proc/cpuinfo") as fh:
            return fh.read()
    except OSError:
        return ""


def run_capture(cmd, env=None, timeout=60):
    """(returncode, stdout, stderr) of a short command. The self-test replaces this."""
    proc = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout)
    return proc.returncode, proc.stdout, proc.stderr


def smi_env(environ):
    """Environment for nvidia-smi and the child alike: PCI bus order, so that
    `--gpu N` names the same board for both (CUDA's default order is fastest-first)."""
    env = dict(environ)
    env["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    return env


def query_gpu(visible, environ):
    """{'name', 'driver', 'cuda'} of the first visible board, from nvidia-smi."""
    first = visible.split(",")[0].strip()
    env = smi_env(environ)
    try:
        rc, out, err = run_capture(["nvidia-smi", "--query-gpu=name,driver_version", "--format=csv,noheader",
                                    "-i", first], env)
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise UsageError(f"cannot run nvidia-smi ({exc}); pass --arch to name the architecture")
    if rc != 0 or not out.strip():
        raise UsageError(f"nvidia-smi -i {first} failed ({err.strip()}); pass --arch to name the architecture")
    fields = [f.strip() for f in out.strip().splitlines()[0].split(",")]
    info = {"name": fields[0], "driver": fields[1] if len(fields) > 1 else "unknown", "cuda": "unknown"}
    try:
        rc, out, _ = run_capture(["nvidia-smi"], env)
        m = re.search(r"CUDA Version:\s*([0-9.]+)", out)
        if rc == 0 and m:
            info["cuda"] = m.group(1)
    except (OSError, subprocess.TimeoutExpired):
        pass
    return info


def machine_identity(info, visible, arch_override, environ):
    """Everything about the machine a golden should record, and the file-name slug."""
    cpuinfo = read_cpuinfo()
    ident = {"host": socket.gethostname(), "cpu_model": cpu_model_name(cpuinfo) or platform.processor() or "unknown",
             "driver": "-", "cuda": "-"}
    if info["device"]:
        try:
            gpu = query_gpu(visible, environ)
        except UsageError:
            if not arch_override:
                raise
            gpu = {"name": "unknown", "driver": "unknown", "cuda": "unknown"}
        ident.update(gpu_name=gpu["name"], driver=gpu["driver"], cuda=gpu["cuda"])
        ident["arch"] = arch_override or gpu_slug(gpu["name"])
    else:
        ident["gpu_name"] = "none (host build)"
        ident["arch"] = arch_override or cpu_slug(cpuinfo)
    return ident


def git_commit():
    """The commit of the checkout this script runs in."""
    try:
        sha = subprocess.run(["git", "-C", REPO, "rev-parse", "--short", "HEAD"], capture_output=True,
                             text=True, timeout=30).stdout.strip()
        dirty = subprocess.run(["git", "-C", REPO, "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, timeout=30).stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        return "unknown"
    if not sha:
        return "unknown"
    return sha + ("-dirty" if dirty else "")


# ----------------------------------------------------------------------------
# running a case
# ----------------------------------------------------------------------------

def flatten(obj, prefix=""):
    """Nested JSON numbers -> {"a.b.c": value}."""
    out = {}
    for key, val in obj.items():
        path = f"{prefix}.{key}" if prefix else key
        if isinstance(val, dict):
            out.update(flatten(val, path))
        elif isinstance(val, (int, float)) and not isinstance(val, bool):
            out[path] = float(val)
    return out


def with_derived(flat):
    out = dict(flat)
    for name, (num, den) in DERIVED.items():
        if num in flat and den in flat and flat[den] > 0:
            out[name] = flat[num] / flat[den]
    return out


def extract_energy(doc):
    """The energy a case is judged on: the SCF total, plus the MP2 correlation
    for --rimp2 (terco's 'total energy' is the SCF's unless a double hybrid)."""
    if "total energy" not in doc:
        return None
    e = float(doc["total energy"])
    if "mp2 opposite-spin energy" in doc and "pt2 opposite-spin fraction" not in doc:
        e += float(doc.get("mp2 opposite-spin energy", 0.0)) + float(doc.get("mp2 same-spin energy", 0.0))
    return e


def case_command(terco, case, json_path):
    return [terco, case["geometry_path"]] + list(case["args"]) + ["--json", json_path, "--silent"]


def child_env(visible, environ):
    env = dict(environ)
    if visible is not None:
        env["CUDA_VISIBLE_DEVICES"] = visible
        env["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
    if not env.get("TERCO_BASIS_PATH"):
        env["TERCO_BASIS_PATH"] = os.path.join(REPO, "basis_sets") + os.pathsep + os.path.join(REPO, "test", "basis_sets")
    return env


def kill_group(proc):
    """SIGKILL the child's whole process group, so no GPU process outlives the bench."""
    try:
        os.killpg(proc.pid, signal.SIGKILL)
    except (ProcessLookupError, PermissionError):
        pass


def run_child(cmd, env, cwd, log_path, timeout):
    """Run cmd in its own session with output in log_path; return its exit status.
    On a timeout (or an interrupt) the whole process group is killed. Raises
    subprocess.TimeoutExpired after the kill."""
    with open(log_path, "w") as log:
        proc = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT, env=env, cwd=cwd, start_new_session=True)
        try:
            return proc.wait(timeout=timeout)
        except BaseException:
            kill_group(proc)
            proc.wait()
            raise


def run_once(terco, case, visible, outdir, tag, max_seconds, environ):
    """One run: returns the parsed JSON. Output is kept in outdir either way.
    Raises RunError (verdict ERROR or NOT CONVERGED) with the log path."""
    stem = os.path.join(outdir, f"{case['name']}.{tag}")
    json_path = stem + ".json"
    log_path = stem + ".log"
    if os.path.exists(json_path):
        os.remove(json_path)    # a run that writes none must not read the previous one's
    cmd = case_command(terco, case, json_path)
    try:
        rc = run_child(cmd, child_env(visible, environ), outdir, log_path, max_seconds)
    except subprocess.TimeoutExpired:
        raise RunError(f"exceeded max_seconds={max_seconds} ({tag}); process group killed. Raise max_seconds in "
                       f"cases.toml only if the case is meant to be that big", log=log_path)
    except OSError as exc:
        raise RunError(f"cannot run {cmd[0]}: {exc}", log=log_path)
    doc = None
    try:
        with open(json_path) as fh:
            doc = json.load(fh)
    except (OSError, ValueError):
        doc = None
    if rc == 2 and doc is not None:
        raise RunError(f"SCF not converged ({tag}): terco exited 2 (JSON kept)", "NOT CONVERGED", log_path)
    if rc != 0:
        raise RunError(f"terco exited with status {rc} ({tag})", log=log_path)
    if doc is None:
        raise RunError(f"no usable JSON at {json_path} ({tag})", log=log_path)
    if doc.get("converged") is False:
        raise RunError(f"SCF reported converged=false ({tag})", "NOT CONVERGED", log_path)
    return doc


def run_case(terco, case, visible, outdir, repeats, warmup, limit, environ, say=print):
    """(docs, None) or (None, RunError). Warm-up runs are run first and discarded."""
    docs = []
    try:
        for w in range(1, warmup + 1):
            say(f"  running {case['name']} (warm-up {w}/{warmup}, discarded) ...")
            run_once(terco, case, visible, outdir, f"warmup{w}", limit, environ)
        for rep in range(1, repeats + 1):
            say(f"  running {case['name']} ({rep}/{repeats}) ...")
            docs.append(run_once(terco, case, visible, outdir, str(rep), limit, environ))
    except RunError as exc:
        return None, exc
    return docs, None


def aggregate(case, docs):
    """Per-repeat JSONs -> {stage: [values]} plus the per-repeat energies,
    iteration counts and gradient RMS values."""
    stages = {}
    for doc in docs:
        flat = with_derived(flatten(doc.get("timings", {})))
        for stage in case["stages"]:
            if stage in flat:
                stages.setdefault(stage, []).append(flat[stage])
    energies = [e for e in (extract_energy(d) for d in docs) if e is not None]
    iters = [int(d["iterations"]) for d in docs if "iterations" in d]
    grads = [float(d["gradient rms"]) for d in docs if "gradient rms" in d]
    return {
        "stages": stages,
        "energies": energies, "iterations_all": iters, "gradients": grads,
        "energy": statistics.median(energies) if energies else None,
        "iterations": iters[0] if iters else None,
        "gradient_rms": statistics.median(grads) if grads else None,
    }


def rows_from_result(case_name, agg, prov, repeats=None):
    rows = []
    for stage, vals in agg["stages"].items():
        row = {
            "case": case_name, "stage": stage,
            "median_s": statistics.median(vals), "min_s": min(vals), "max_s": max(vals),
            "repeats": len(vals),
            "energy_hartree": agg["energy"], "iterations": agg["iterations"],
            "gradient_rms": agg["gradient_rms"],
        }
        row.update(prov)
        rows.append(row)
    return rows


# ----------------------------------------------------------------------------
# goldens
# ----------------------------------------------------------------------------

def _fmt(key, val):
    if val is None:
        return ""
    if key in ("median_s", "min_s", "max_s"):
        return f"{val:.6f}"    # microseconds: the DF stages of a water hexamer are milliseconds
    if key == "energy_hartree":
        return f"{val:.10f}"
    if key == "gradient_rms":
        return f"{val:.10e}"   # compared at gradient_tol (1e-6 by default): three digits were not enough
    return str(val)


def sort_rows(rows, case_order):
    """cases.toml order, then stage name; cases cases.toml does not know follow, by name."""
    return sorted(rows, key=lambda r: (case_order.get(r["case"], len(case_order)), r["case"], r["stage"]))


def write_golden(path, rows, case_order=None):
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    rows = sort_rows(rows, case_order or {})
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=GOLDEN_FIELDS, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({k: _fmt(k, row.get(k)) for k in GOLDEN_FIELDS})


def read_golden(path):
    """Golden CSV -> list of row dicts with numbers converted; [] if absent."""
    if not os.path.isfile(path):
        return []
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh)
        missing = [f for f in GOLDEN_FIELDS if f not in (reader.fieldnames or [])]
        if missing:
            raise UsageError(f"{path}: missing columns {missing}; it was written by another version of run.py "
                             f"and has to be recorded again")
        for raw in reader:
            row = dict(raw)
            for key in ("median_s", "min_s", "max_s", "energy_hartree", "gradient_rms"):
                row[key] = float(raw[key]) if raw[key] != "" else None
            for key in ("repeats", "iterations"):
                row[key] = int(raw[key]) if raw[key] != "" else None
            rows.append(row)
    return rows


def merge_golden(old_rows, new_rows):
    """The new rows replace every old row of the cases that ran; other cases stay
    (with their own provenance)."""
    new_cases = {r["case"] for r in new_rows}
    kept = [r for r in old_rows if r["case"] not in new_cases]
    return kept + new_rows


# ----------------------------------------------------------------------------
# judging
# ----------------------------------------------------------------------------

def judge_stage(golden_med, now_med, fail, warn, min_seconds, judged=True):
    """(ratio, verdict) for one stage.
    INFO: not a judged stage; reported, never fails.
    SMALL: both the golden and the new time are under min_seconds (a 40x growth
    from 0.04 s to 1.6 s is judged, not hidden here).
    FAST: more than `fail` times faster (information only)."""
    if golden_med is None or now_med is None:
        return None, "MISSING"
    ratio = (now_med / golden_med) if golden_med > 0 else None
    if not judged:
        return ratio, "INFO"
    if max(golden_med, now_med) < min_seconds:
        return ratio, "SMALL"
    if ratio is None:
        return None, "FAIL"    # from nothing to over min_seconds
    if ratio > fail:
        return ratio, "FAIL"
    if ratio > warn:
        return ratio, "WARN"
    if ratio < 1.0 / fail:
        return ratio, "FAST"
    return ratio, "OK"


def _num(val, digits=3):
    return "-" if val is None else f"{val:.{digits}f}"


def _spread(vals_min, vals_max):
    return f"{vals_min:.3f}-{vals_max:.3f}"


def consistency_rows(name, case, agg, defaults):
    """Rows for repeats that disagree with each other: energy, gradient RMS, SCF iterations."""
    rows = []
    etol = case.get("energy_tol", defaults["energy_tol"])
    es = agg["energies"]
    if len(es) > 1 and max(es) - min(es) > etol:
        rows.append((name, "energy [Eh] repeats", "-", f"{min(es):.10f}", f"d={max(es) - min(es):.1e}", "-", "FAIL",
                     f"repeats disagree, tol {etol:.0e}"))
    gtol = case.get("gradient_tol", defaults["gradient_tol"])
    gs = agg["gradients"]
    if len(gs) > 1 and max(gs) - min(gs) > gtol:
        rows.append((name, "gradient rms repeats", "-", f"{min(gs):.8f}", f"d={max(gs) - min(gs):.1e}", "-", "FAIL",
                     f"repeats disagree, tol {gtol:.0e}"))
    its = agg["iterations_all"]
    if len(set(its)) > 1:
        # A converged SCF near its threshold takes n or n + 1 iterations from one run to the next on a
        # GPU (measured: gly100 6-31G, 15 or 16). The per-build means do not depend on the count, so this
        # is the case's iteration verdict (WARN unless the case says FAIL), not a failure of its own.
        rows.append((name, "scf iterations repeats", "-", f"{min(its)}-{max(its)}", "-", "-", case["iteration_verdict"],
                     "repeats took different iteration counts"))
    return rows


def compare(cases_run, results, golden_rows, defaults, per_stage, have_golden, all_case_names=None):
    """(table, hints). Table rows are (case, stage, golden, now, ratio, warn/fail, verdict, note)
    for the cases that ran. `results[name]` is {"agg":..., "rows": [...]}."""
    table, hints = [], set()
    by_key = {(r["case"], r["stage"]): r for r in golden_rows}
    golden_cases = {r["case"] for r in golden_rows}
    for case in cases_run:
        name = case["name"]
        now_rows = {r["stage"]: r for r in results[name]["rows"]}
        agg = results[name]["agg"]
        table.extend(consistency_rows(name, case, agg, defaults))
        if not have_golden or name not in golden_cases:
            for stage in case["stages"]:
                r = now_rows.get(stage)
                if r is None:
                    table.append((name, stage, "-", "-", "-", "-", "WARN", "missing now"))
                else:
                    table.append((name, stage, "-", _num(r["median_s"]), "-", "-",
                                  "no golden" if not have_golden else "NEW", _spread(r["min_s"], r["max_s"])))
            continue
        g_any = next(r for r in golden_rows if r["case"] == name)
        # energy: every repeat against the golden, so a faster wrong answer fails
        etol = case.get("energy_tol", defaults["energy_tol"])
        es, e_gold = agg["energies"], g_any["energy_hartree"]
        if not es or e_gold is None:
            table.append((name, "energy [Eh]", _num(e_gold, 10), _num(None), "-", "-", "WARN", "energy missing"))
        else:
            worst = max(es, key=lambda e: abs(e - e_gold))
            bad = sum(1 for e in es if abs(e - e_gold) > etol)
            table.append((name, "energy [Eh]", f"{e_gold:.10f}", f"{worst:.10f}", f"d={abs(worst - e_gold):.1e}", "-",
                          "FAIL" if bad else "OK", f"tol {etol:.0e}" + (f", {bad}/{len(es)} repeats beyond" if bad else "")))
        gs, g_gold = agg["gradients"], g_any["gradient_rms"]
        if gs and g_gold is not None:
            gtol = case.get("gradient_tol", defaults["gradient_tol"])
            worst = max(gs, key=lambda g: abs(g - g_gold))
            bad = sum(1 for g in gs if abs(g - g_gold) > gtol)
            table.append((name, "gradient rms", f"{g_gold:.8f}", f"{worst:.8f}", f"d={abs(worst - g_gold):.1e}", "-",
                          "FAIL" if bad else "OK", f"tol {gtol:.0e}" + (f", {bad}/{len(gs)} repeats beyond" if bad else "")))
        # iteration count: its own row, never folded into a timing
        its, it_gold = agg["iterations_all"], g_any["iterations"]
        if its and it_gold is not None:
            differ = any(i != it_gold for i in its)
            table.append((name, "scf iterations", str(it_gold), str(its[0]), "-", "-",
                          case["iteration_verdict"] if differ else "OK", ""))
        judged = set(case["judge"])
        for stage in case["stages"]:
            thr = thresholds_for(defaults, per_stage, case, stage)
            g, n = by_key.get((name, stage)), now_rows.get(stage)
            if g is None and n is None:
                continue
            if g is None:
                table.append((name, stage, "-", _num(n["median_s"]), "-", "-", "NEW", _spread(n["min_s"], n["max_s"])))
                continue
            if n is None:
                table.append((name, stage, _num(g["median_s"]), "-", "-", "-", "WARN", "in the golden, missing now"))
                continue
            fail, warn, source = effective_ratios(thr, g)
            ratio, verdict = judge_stage(g["median_s"], n["median_s"], fail, warn, thr["min_seconds"], stage in judged)
            if verdict == "FAST":
                hints.add("FAST")
            table.append((name, stage, _num(g["median_s"]), _num(n["median_s"]),
                          "-" if ratio is None else f"{ratio:.3f}",
                          f"{warn:.2f}/{fail:.2f}" + ("*" if source == "spread" else "") if stage in judged else "-",
                          verdict, _spread(n["min_s"], n["max_s"])))
        for stage in sorted(set(r["stage"] for r in golden_rows if r["case"] == name) - set(case["stages"])):
            g = by_key[(name, stage)]
            table.append((name, stage, _num(g["median_s"]), "-", "-", "-", "WARN", "in the golden, no longer in cases.toml"))
    if all_case_names is not None:
        for gone in sorted(golden_cases - set(all_case_names)):
            table.append((gone, "(case)", "-", "-", "-", "-", "WARN", "in the golden, no longer in cases.toml"))
    return table, hints


def exit_status_of(table):
    if any(row[6] in ERROR_VERDICTS for row in table):
        return 2
    return 1 if any(row[6] == "FAIL" for row in table) else 0


HEADER = ("case", "stage", "golden s", "now s", "ratio", "warn/fail", "verdict", "now min-max s / note")


def format_table(table):
    rows = [HEADER] + [tuple(str(c) for c in r) for r in table]
    widths = [max(len(r[i]) for r in rows) for i in range(len(HEADER))]
    lines = []
    for k, r in enumerate(rows):
        lines.append("  ".join(c.ljust(widths[i]) if i in (0, 1, 6) else c.rjust(widths[i]) if i < 6 else c
                               for i, c in enumerate(r)).rstrip())
        if k == 0:
            lines.append("  ".join("-" * w for w in widths))
    return "\n".join(lines)


def summary_line(table, n_cases, skipped):
    counts = {}
    for row in table:
        counts[row[6]] = counts.get(row[6], 0) + 1
    parts = ", ".join(f"{v} {k}" for k, v in sorted(counts.items()))
    code = exit_status_of(table)
    status = {0: "no regression", 1: "REGRESSION", 2: "INCOMPLETE (a case errored, did not converge or was refused)"}[code]
    return f"summary: {status}: {n_cases} case(s) run, {len(skipped)} skipped; {parts or 'no rows'}"


FAST_HINT = ("hint: a judged stage is more than its fail ratio FASTER than the golden. If that is an intended "
             "speed-up, re-record the golden in a PR that states it; otherwise a later regression hides under the "
             "stale golden.")


def describe_changes(old_rows, new_rows):
    """What --record changes versus the previous golden, for printing."""
    old = {(r["case"], r["stage"]): r for r in old_rows}
    lines = []
    seen_energy = set()
    for r in new_rows:
        key = (r["case"], r["stage"])
        if key not in old:
            lines.append(f"  new      {r['case']} {r['stage']}: {r['median_s']:.3f} s")
        else:
            o = old[key]
            if o["median_s"] and o["median_s"] > 0:
                lines.append(f"  changed  {r['case']} {r['stage']}: {o['median_s']:.3f} -> {r['median_s']:.3f} s "
                             f"(x{r['median_s'] / o['median_s']:.3f})")
            if (r["case"] not in seen_energy and o["energy_hartree"] is not None and r["energy_hartree"] is not None
                    and abs(o["energy_hartree"] - r["energy_hartree"]) > 1e-10):
                seen_energy.add(r["case"])
                lines.append(f"  ENERGY   {r['case']}: {o['energy_hartree']:.10f} -> {r['energy_hartree']:.10f} Eh "
                             f"(d={abs(o['energy_hartree'] - r['energy_hartree']):.1e})")
    new_keys = {(r["case"], r["stage"]) for r in new_rows}
    new_cases = {r["case"] for r in new_rows}
    for key, o in old.items():
        if key[0] in new_cases and key not in new_keys:
            lines.append(f"  removed  {key[0]} {key[1]} (was {o['median_s']:.3f} s)")
    return lines


def machine_warnings(golden_rows, cases, ident):
    """Same family, different board or host: say so (judged device stages are still compared;
    host-side INFO stages will differ)."""
    warnings = []
    names = {c["name"] for c in cases}
    seen = set()
    for r in golden_rows:
        if r["case"] not in names:
            continue
        for key, now, label in (("gpu_name", ident.get("gpu_name"), "GPU"), ("cpu_model", ident.get("cpu_model"), "CPU")):
            then = r.get(key)
            if then and now and then != now and (key, then) not in seen:
                seen.add((key, then))
                warnings.append(f"the golden was recorded on {label} '{then}' (host {r.get('host')}), this machine has "
                                f"'{now}': judged device stages are still compared, host-side INFO stages will differ")
    return warnings


# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------

def build_parser():
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter,
                                epilog="Exit status: 0 no regression, 1 a FAIL, 2 usage error or an errored case.")
    p.add_argument("--terco", help="the terco binary (default: the newest build*/app/terco in the repo)")
    p.add_argument("--cases", help="comma-separated case names (default: all the build can run)")
    p.add_argument("--quick", action="store_true", help="only the cases marked quick in cases.toml")
    p.add_argument("--repeats", type=int, default=3, help="runs per case; the median is judged (default 3; "
                                                           "record with 5 or more for a spread-derived threshold)")
    p.add_argument("--warmup", type=int, default=None,
                   help="extra runs per case, discarded, before the measured ones (default 1 for a GPU build, 0 for a host build)")
    p.add_argument("--gpu", type=int, help="physical GPU index in PCI-bus order (sets CUDA_VISIBLE_DEVICES and "
                                           "CUDA_DEVICE_ORDER=PCI_BUS_ID); required for a GPU build unless "
                                           "CUDA_VISIBLE_DEVICES is already set")
    p.add_argument("--allow-gpu0", action="store_true", help="permit GPU 0 (refused otherwise)")
    p.add_argument("--arch", help="golden file slug (default: the GPU family from nvidia-smi, or cpu-<model>)")
    p.add_argument("--goldens", help="directory of golden CSVs (default: bench/regression/goldens)")
    p.add_argument("--config", help="cases file (default: bench/regression/cases.toml)")
    p.add_argument("--out", help="save the raw JSON and log of every run, and result.csv, in this directory")
    p.add_argument("--record", action="store_true", help="write/update this architecture's golden from this run")
    p.add_argument("--allow-build-mismatch", action="store_true",
                   help="compare or record against a golden whose build differs (device, GPU arch, LMAX, DF LMAX, "
                        "compiler family, ...); refused otherwise")
    p.add_argument("--list", action="store_true", help="list the cases and exit")
    p.add_argument("--dry-run", action="store_true", help="print the commands and exit")
    p.add_argument("--self-test", action="store_true", help="run the script's own unit tests (no terco needed)")
    return p


def select_cases(all_cases, args):
    if args.cases and args.quick:
        raise UsageError("--cases and --quick are exclusive")
    if args.cases:
        names = [n.strip() for n in args.cases.split(",") if n.strip()]
        known = {c["name"] for c in all_cases}
        bad = [n for n in names if n not in known]
        if bad:
            raise UsageError(f"unknown case(s) {', '.join(bad)}; --list shows them")
        return [c for c in all_cases if c["name"] in names]
    if args.quick:
        return [c for c in all_cases if c["quick"]]
    return list(all_cases)


def list_cases(cases):
    for c in cases:
        flags = []
        if c["quick"]:
            flags.append("quick")
        if c["needs"]:
            flags.append("needs " + "+".join(c["needs"]))
        print(f"{c['name']:<14} {c['geometry']:<26} {' '.join(c['args'])}"
              f"{'   [' + ', '.join(flags) + ']' if flags else ''}")
        print(f"{'':<14} {c['description']}")
        print(f"{'':<14} judged: {', '.join(c['judge']) or '(none)'}")


def main(argv=None, environ=None):
    environ = os.environ if environ is None else environ
    args = build_parser().parse_args(argv)
    if args.self_test:
        return self_test()
    if args.repeats < 1:
        raise UsageError("--repeats must be at least 1")
    if args.warmup is not None and args.warmup < 0:
        raise UsageError("--warmup cannot be negative")
    if args.gpu is not None and args.gpu < 0:
        raise UsageError(f"--gpu {args.gpu}: a GPU index cannot be negative")
    # The child runs in the output directory, so every path it or we use must not depend on the caller's cwd.
    config_path = os.path.abspath(args.config or os.path.join(HERE, "cases.toml"))
    goldens_dir = os.path.abspath(args.goldens or os.path.join(HERE, "goldens"))
    out_arg = os.path.abspath(args.out) if args.out else None
    terco_arg = os.path.abspath(args.terco) if args.terco else None
    defaults, per_stage, all_cases = load_config(config_path)
    selected = select_cases(all_cases, args)
    if args.list:
        list_cases(selected)
        return 0
    terco = discover_terco(terco_arg)
    info = probe_build(terco)
    caps = build_capabilities(info)
    warmup = args.warmup if args.warmup is not None else (1 if info["device"] else 0)
    cases_run, skipped = [], []
    for c in selected:
        missing = [n for n in c["needs"] if n not in caps]
        if missing:
            skipped.append((c["name"], missing))
            print(f"note: skipping {c['name']}: needs {', '.join(missing)} "
                  f"(this build: device={'yes' if info['device'] else 'no'}, lmax={info.get('lmax')})")
        else:
            cases_run.append(c)
    if not cases_run:
        raise UsageError("no case can run with this build")
    visible = resolve_gpu(info["device"], args.gpu, args.allow_gpu0, environ)
    if args.dry_run:
        env = f"CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES={visible} " if visible is not None else ""
        for c in cases_run:
            cmd = case_command(terco, c, "<out>/" + c["name"] + ".<run>.json")
            print(f"{env}{' '.join(cmd)}   # {warmup} warm-up + {args.repeats} repeat(s), "
                  f"timeout {c.get('max_seconds', defaults['max_seconds'])} s")
        return 0
    ident = machine_identity(info, visible, args.arch, environ)
    arch = ident["arch"]
    golden_path = os.path.join(goldens_dir, arch + ".csv")
    outdir = out_arg or tempfile.mkdtemp(prefix="terco-bench-")
    os.makedirs(outdir, exist_ok=True)
    build_str = build_string(info)
    omp, mkl = environ.get("OMP_NUM_THREADS"), environ.get("MKL_NUM_THREADS")
    prov = {
        "binary_commit": info.get("commit", "unknown"), "checkout_commit": git_commit(),
        "date": datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "host": ident["host"], "cpu_model": ident["cpu_model"], "gpu_name": ident["gpu_name"],
        "driver": ident["driver"], "cuda": ident["cuda"],
        "omp_threads": omp if omp else "unset", "mkl_threads": mkl if mkl else "unset",
        "arch": arch, "build": build_str, "compiler": info.get("compiler", "unknown"),
    }
    print(f"terco {terco} (commit {prov['binary_commit']}, from the last configure)\n"
          f"arch {arch} ({ident['gpu_name']}; driver {ident['driver']}, CUDA {ident['cuda']}); cpu {ident['cpu_model']}\n"
          f"build {build_str}; compiler {prov['compiler']}\n"
          f"{warmup} warm-up + {args.repeats} repeat(s); OMP_NUM_THREADS={prov['omp_threads']} "
          f"MKL_NUM_THREADS={prov['mkl_threads']}; raw output in {outdir}")
    warnings = []
    if not omp:
        print("note: OMP_NUM_THREADS is unset, so host-side stages use the library's default thread count; "
              "set it for comparable host-side numbers")
    if gpu_is_uuid(visible):
        warnings.append(f"CUDA_VISIBLE_DEVICES={visible} is a GPU UUID: recorded as given, but it cannot be checked "
                        f"against the 'never GPU 0' rule")
    golden_rows = read_golden(golden_path)
    # A golden is only meaningful against the build it was recorded with: refuse per case, before running.
    refused, table_pre = {}, []
    for c in cases_run:
        grows = [r for r in golden_rows if r["case"] == c["name"]]
        if not grows:
            continue
        diff = build_mismatch(grows[0].get("build"), grows[0].get("compiler"), build_str, prov["compiler"])
        if diff:
            if args.allow_build_mismatch:
                warnings.append(f"case {c['name']}: build differs from the golden ({'; '.join(diff)}); compared "
                                f"anyway (--allow-build-mismatch)")
            else:
                refused[c["name"]] = diff
                table_pre.append((c["name"], "(build)", "-", "-", "-", "-", "BUILD MISMATCH",
                                  "; ".join(diff) + ". Rebuild like the golden, record a new golden, or pass "
                                                    "--allow-build-mismatch"))
    runnable = [c for c in cases_run if c["name"] not in refused]
    results, errors = {}, []
    for c in runnable:
        limit = c.get("max_seconds", defaults["max_seconds"])
        docs, err = run_case(terco, c, visible, outdir, args.repeats, warmup, limit, environ)
        if err is not None:
            print(f"  ERROR in {c['name']}: {err}", file=sys.stderr)
            errors.append((c["name"], err))
            continue
        agg = aggregate(c, docs)
        results[c["name"]] = {"agg": agg, "rows": rows_from_result(c["name"], agg, prov)}
    ok_cases = [c for c in runnable if c["name"] in results]
    all_names = {c["name"] for c in all_cases} if not (args.cases or args.quick) else None
    table, hints = compare(ok_cases, results, golden_rows, defaults, per_stage, bool(golden_rows), all_names)
    table = table_pre + table
    for name, err in errors:
        note = str(err) + (f" [log {err.log}]" if err.log else "")
        table.append((name, "(run)", "-", "-", "-", "-", err.verdict, note))
    case_order = {c["name"]: i for i, c in enumerate(all_cases)}
    new_rows = [r for c in ok_cases for r in results[c["name"]]["rows"]]
    write_golden(os.path.join(outdir, "result.csv"), new_rows, case_order)
    warnings.extend(machine_warnings(golden_rows, ok_cases, ident))
    print()
    print(format_table(table))
    print()
    print(summary_line(table, len(ok_cases), skipped))
    if "FAST" in hints:
        print(FAST_HINT)
    for w in warnings:
        print(f"warning: {w}")
    if not golden_rows:
        print(f"no golden for {arch} ({golden_path}): run with --record to create it from this run")
    if args.record:
        if args.repeats < SPREAD_REPEATS:
            print(f"note: recording with {args.repeats} repeat(s); a spread-derived threshold needs at least "
                  f"{SPREAD_REPEATS}, so this golden is compared with the configured ratios only")
        changes = describe_changes(golden_rows, new_rows)
        print(f"\n--record: {golden_path}")
        print("\n".join(changes) if changes else "  (nothing to compare with: first record)")
        skipped_rec = [n for n, _ in errors] + sorted(refused)
        if skipped_rec:
            print(f"  not recorded (errored or refused; their old rows stay): {', '.join(skipped_rec)}")
        if new_rows:
            write_golden(golden_path, merge_golden(golden_rows, new_rows), case_order)
            print(f"wrote {golden_path} ({len(new_rows)} rows from {len(ok_cases)} case(s)). "
                  f"Commit it only for an intended change, and say so in the PR.")
    return exit_status_of(table)


# ----------------------------------------------------------------------------
# self test
# ----------------------------------------------------------------------------

FAKE_TERCO = r'''#!{python}
import json, os, sys, subprocess
args = sys.argv[1:]
if "--version" in args:
    print("terco version 0.0.0")
    print("build: device={device} {gpu_arch}lmax={lmax} df_lmax=2 df_lmax_aux=4 cart_order=lexicographic mpi=no "
          "multicore=no commit=abc1234 compiler={compiler}")
    sys.exit(0)
counter = os.environ.get("FAKE_COUNTER")
call = 1
if counter:
    call = int(open(counter).read()) + 1 if os.path.exists(counter) else 1
    open(counter, "w").write(str(call))
if os.environ.get("FAKE_FAIL_ARG") and os.environ["FAKE_FAIL_ARG"] in args:
    print("boom"); sys.exit(3)
if os.environ.get("FAKE_HANG"):
    if os.environ.get("FAKE_SPAWN"):
        child = subprocess.Popen(["sleep", "60"])
        open(os.environ["FAKE_SPAWN"], "w").write(str(child.pid))
    import time; time.sleep(60)
scale = float(os.environ.get("FAKE_SCALE", "1.0"))
if os.environ.get("FAKE_SLOWFIRST") and call == 1:
    scale *= 5.0
iters = int(os.environ.get("FAKE_ITERS", "10"))
if os.environ.get("FAKE_ITERS_AT") and call == int(os.environ["FAKE_ITERS_AT"]):
    iters += 1
shift = float(os.environ.get("FAKE_ESHIFT", "0.0"))
if os.environ.get("FAKE_ESHIFT_AT") and call == int(os.environ["FAKE_ESHIFT_AT"]):
    shift += float(os.environ.get("FAKE_ESHIFT_ONE", "1e-5"))
# per-call jitter, deterministic: the 5 recorded repeats get a known spread
jit = float(os.environ.get("FAKE_JITTER", "0.0")) * ((call % 5) - 2)
doc = {{"total energy": -76.0 + shift, "iterations": iters, "converged": True, "gradient rms": 0.01,
       "timings": {{"total_wall": 2.0 * scale, "scf": {{"wall": 1.5 * scale * (1 + jit), "fock_mean_after_first": 0.1 * scale * (1 + jit),
                                                    "fock_total": 1.0 * scale, "n_fock_builds": iters, "n_iterations": iters}},
                   "gradient": {{"total": 0.2 * scale}}}}}}
if os.environ.get("FAKE_NOCONV"):
    doc["converged"] = False
open(args[args.index("--json") + 1], "w").write(json.dumps(doc))
sys.exit(2 if os.environ.get("FAKE_NOCONV") else 0)
'''

FAKE_CASES = """
[defaults]
max_seconds = 3
[[case]]
name = "a"
description = "synthetic"
geometry = "geom.xyz"
args = ["-b", "x"]
stages = ["scf.wall", "scf.fock_mean_after_first", "scf.fock_per_iteration", "gradient.total", "total_wall"]
judge = ["scf.fock_mean_after_first", "gradient.total"]
quick = true
[[case]]
name = "b"
description = "synthetic, needs f"
geometry = "geom.xyz"
args = ["-b", "y"]
needs = ["lmax3"]
stages = ["scf.wall"]
judge = ["scf.wall"]
[[case]]
name = "c"
description = "synthetic, fails when FAKE_FAIL_ARG=z"
geometry = "geom.xyz"
args = ["-b", "z"]
stages = ["scf.wall", "scf.fock_mean_after_first"]
judge = ["scf.fock_mean_after_first"]
"""


class _Checker:
    def __init__(self):
        self.failures = []
        self.count = 0

    def check(self, cond, msg):
        self.count += 1
        if not cond:
            self.failures.append(msg)
            print(f"  FAIL: {msg}")


def self_test():
    t = _Checker()
    check = t.check

    # flatten / aggregate / energy
    flat = flatten({"total_wall": 2, "scf": {"wall": 1.5, "n_iterations": 4}, "x": True, "s": "str"})
    check(flat == {"total_wall": 2.0, "scf.wall": 1.5, "scf.n_iterations": 4.0}, f"flatten {flat}")
    check(with_derived({"scf.fock_total": 3.0, "scf.n_fock_builds": 6.0})["scf.fock_per_iteration"] == 0.5, "derived stage")
    check("scf.fock_per_iteration" not in with_derived({"scf.fock_total": 3.0}), "derived stage needs its inputs")
    check(extract_energy({"total energy": -1.0}) == -1.0, "energy plain")
    e = extract_energy({"total energy": -1.0, "mp2 opposite-spin energy": -0.2, "mp2 same-spin energy": -0.1})
    check(abs(e + 1.3) < 1e-14, f"energy with mp2 {e}")
    check(extract_energy({"total energy": -1.0, "mp2 opposite-spin energy": -0.2, "mp2 same-spin energy": -0.1,
                          "pt2 opposite-spin fraction": 0.5}) == -1.0, "double hybrid total is already complete")
    case = {"name": "c", "stages": ["scf.wall", "nothere"]}
    docs = [{"total energy": -1.0, "iterations": 7, "timings": {"scf": {"wall": v}}} for v in (3.0, 1.0, 5.0)]
    agg = aggregate(case, docs)
    check(agg["stages"] == {"scf.wall": [3.0, 1.0, 5.0]}, f"aggregate {agg['stages']}")
    prov0 = {k: "p" for k in PROVENANCE}
    rows = rows_from_result("c", agg, prov0)
    check(len(rows) == 1 and rows[0]["median_s"] == 3.0 and rows[0]["min_s"] == 1.0 and rows[0]["max_s"] == 5.0,
          "median (3.0), not the mean (3.0 here)")
    skewed = aggregate(case, [{"timings": {"scf": {"wall": v}}} for v in (1.0, 1.0, 10.0)])
    check(rows_from_result("c", skewed, prov0)[0]["median_s"] == 1.0, "median versus mean (mean would be 4.0)")

    # verdicts, exactly at the ratios (10 -> 10.5 is exactly 1.05: not above warn)
    F, W, M = 1.10, 1.05, 0.05
    check(judge_stage(10.0, 10.4, F, W, M)[1] == "OK", "ratio 1.04 OK")
    check(judge_stage(100.0, 105.0, F, W, M)[1] == "OK", "exactly the warn ratio is OK")
    check(judge_stage(100.0, 105.1, F, W, M)[1] == "WARN", "just above warn")
    check(judge_stage(100.0, 110.0, F, W, M)[1] == "WARN", "exactly the fail ratio is WARN")
    check(judge_stage(100.0, 110.1, F, W, M)[1] == "FAIL", "just above fail")
    check(judge_stage(10.0, 8.0, F, W, M)[1] == "FAST", "much faster is FAST, not a failure")
    check(judge_stage(0.01, 0.02, F, W, M)[1] == "SMALL", "both below min_seconds is SMALL")
    check(judge_stage(0.04, 2.0, F, W, M)[1] == "FAIL", "0.04 -> 2.0 s (50x) is judged and fails")
    check(judge_stage(0.04, 1.6, F, W, M)[1] == "FAIL", "0.04 -> 1.6 s (40x) is judged and fails")
    check(judge_stage(2.0, 0.04, F, W, M)[1] == "FAST", "shrinking below min_seconds is judged too")
    check(judge_stage(10.0, 100.0, F, W, M, judged=False)[1] == "INFO", "an unjudged stage is INFO, never FAIL")
    check(judge_stage(None, 1.0, F, W, M)[1] == "MISSING", "missing golden")
    check(judge_stage(0.0, 1.0, F, W, M)[1] == "FAIL", "from 0 to 1 s fails")
    defaults = dict(DEFAULTS)
    per_stage = {"scf.wall": {"fail": 1.5}}
    c = {"thresholds": {"scf.wall": {"warn": 1.2}}}
    tt = thresholds_for(defaults, per_stage, c, "scf.wall")
    check(tt["fail"] == 1.5 and tt["warn"] == 1.2 and tt["min_seconds"] == 0.05, f"threshold layering {tt}")
    tc = thresholds_for({"fail": 1.10, "warn": 1.05, "min_seconds": 0.05}, {}, {"min_seconds": 0.002}, "any.stage")
    check(tc["min_seconds"] == 0.002 and tc["fail"] == 1.10, f"a case-level min_seconds applies to all its stages: {tc}")
    check(abs(float(_fmt("gradient_rms", 0.01912233)) - 0.01912233) < 1e-9, "the gradient RMS keeps its digits in the golden")
    check(_fmt("median_s", 0.0049) == "0.004900", "times are stored to the microsecond")
    check(thresholds_for(defaults, per_stage, {}, "other")["fail"] == 1.10, "global default")

    # spread-derived thresholds
    thr = {"fail": 1.10, "warn": 1.05, "min_seconds": 0.05}
    few = {"repeats": 3, "median_s": 10.0, "min_s": 8.0, "max_s": 12.0}
    check(effective_ratios(thr, few) == (1.10, 1.05, "configured"), "fewer than 5 repeats: configured ratios")
    tight = {"repeats": 5, "median_s": 10.0, "min_s": 9.99, "max_s": 10.01}
    f, w, src = effective_ratios(thr, tight)
    check(f == 1.10 and abs(w - 1.05) < 1e-12 and src == "spread", f"tiny spread keeps fail 1.10, warn halfway: {(f, w, src)}")
    wide = {"repeats": 5, "median_s": 10.0, "min_s": 9.0, "max_s": 11.0}   # spread 0.2 -> 1 + 0.8
    f, w, src = effective_ratios(thr, wide)
    check(abs(f - 1.8) < 1e-12 and abs(w - 1.4) < 1e-12 and src == "spread", f"wide spread: fail 1.8, warn 1.4: {(f, w)}")
    check(effective_ratios(thr, None) == (1.10, 1.05, "configured"), "no golden row: configured")

    # compare: energy, iterations, missing/new, no golden, INFO, thresholds
    prov = {k: "p" for k in PROVENANCE}

    def golden(med, e=-76.0, it=10, repeats=3, name="a", spread=0.0):
        return [dict({"case": name, "stage": s, "median_s": m, "min_s": m * (1 - spread), "max_s": m * (1 + spread),
                      "repeats": repeats, "energy_hartree": e, "iterations": it, "gradient_rms": None}, **prov)
                for s, m in med.items()]

    case_a = {"name": "a", "stages": ["s1", "s2"], "judge": ["s1", "s2"], "iteration_verdict": "WARN"}

    def now(med, e=(-76.0,), it=10, name="a"):
        n = max(len(e), 1)
        agg = {"stages": {s: [m] * n for s, m in med.items()}, "energies": list(e), "iterations_all": [it] * n,
               "gradients": [], "energy": statistics.median(e) if e else None, "iterations": it, "gradient_rms": None}
        return {name: {"agg": agg, "rows": rows_from_result(name, agg, prov)}}

    def cmp(case, res, gold, have=True, names=None):
        return compare([case], res, gold, defaults, {}, have, names)

    g = golden({"s1": 10.0, "s2": 10.0})
    tab, _ = cmp(case_a, now({"s1": 10.1, "s2": 10.2}), g)
    check(exit_status_of(tab) == 0 and all(r[6] == "OK" for r in tab), f"all OK: {tab}")
    tab, _ = cmp(case_a, now({"s1": 10.1, "s2": 12.0}), g)
    check(exit_status_of(tab) == 1, "a slow judged stage fails")
    tab, _ = cmp(dict(case_a, judge=["s1"]), now({"s1": 10.1, "s2": 12.0}), g)
    check(exit_status_of(tab) == 0 and any(r[1] == "s2" and r[6] == "INFO" for r in tab), "an INFO stage never fails")
    tab, _ = cmp(dict(case_a, judge=[]), now({"s1": 30.0, "s2": 30.0}), g)
    check(exit_status_of(tab) == 0, "nothing judged: nothing fails on time")
    tab, _ = cmp(case_a, now({"s1": 5.0, "s2": 5.0}, e=(-76.0 + 1e-6,)), g)
    check(exit_status_of(tab) == 1 and any(r[1].startswith("energy") and r[6] == "FAIL" for r in tab),
          "a faster wrong answer fails")
    tab, _ = cmp(case_a, now({"s1": 10.0, "s2": 10.0}, e=(-76.0, -76.0, -76.0 + 5e-9)), g)
    check(exit_status_of(tab) == 0, "repeats within tolerance pass")
    tab, _ = cmp(case_a, now({"s1": 10.0, "s2": 10.0}, e=(-76.0, -76.0, -76.0 + 1e-6)), g)
    check(exit_status_of(tab) == 1 and any(r[1] == "energy [Eh]" and r[6] == "FAIL" for r in tab),
          "one repeat's wrong energy fails even though the median is right")
    tab, _ = cmp(case_a, now({"s1": 10.0, "s2": 10.0}, e=(-76.0, -76.0 + 1e-6)), g)
    check(any(r[1] == "energy [Eh] repeats" and r[6] == "FAIL" for r in tab), "repeats disagreeing with each other fail")
    tab, _ = cmp(case_a, now({"s1": 10.0, "s2": 10.0}, e=(-76.0 + 1e-10,)), g)
    check(exit_status_of(tab) == 0, "an energy within tolerance passes")
    tab, _ = cmp(dict(case_a, energy_tol=1e-3), now({"s1": 10.0, "s2": 10.0}, e=(-76.0 + 1e-5,)), g)
    check(exit_status_of(tab) == 0, "per-case energy_tol")
    tab, _ = cmp(case_a, now({"s1": 10.0, "s2": 10.0}, it=12), g)
    check(any(r[1] == "scf iterations" and r[6] == "WARN" for r in tab) and exit_status_of(tab) == 0,
          "a changed iteration count is a WARN row")
    tab, _ = cmp(dict(case_a, iteration_verdict="FAIL"), now({"s1": 10.0, "s2": 10.0}, it=12), g)
    check(exit_status_of(tab) == 1, "iteration_verdict = FAIL makes it a failure")
    res = now({"s1": 10.0, "s2": 10.0})
    res["a"]["agg"]["iterations_all"] = [10, 11]
    tab, _ = cmp(case_a, res, g)
    check(any(r[1] == "scf iterations repeats" and r[6] == "WARN" for r in tab), "repeats with different iteration counts warn")
    tab, _ = cmp(dict(case_a, iteration_verdict="FAIL"), res, g)
    check(any(r[1] == "scf iterations repeats" and r[6] == "FAIL" for r in tab),
          "repeats with different iteration counts fail a case whose iteration verdict is FAIL")
    tab, _ = cmp(case_a, now({"s1": 10.0}), g)
    check(any(r[1] == "s2" and r[6] == "WARN" for r in tab) and exit_status_of(tab) == 0, "stage in golden, missing now: WARN")
    tab, _ = cmp(case_a, now({"s1": 10.0, "s2": 10.0, "s3": 1.0}), golden({"s1": 10.0}))
    check(any(r[1] == "s2" and r[6] == "NEW" for r in tab), "stage now, not in golden: NEW")
    tab, _ = cmp(case_a, now({"s1": 10.0, "s2": 10.0}), [], have=False)
    check(all(r[6] == "no golden" for r in tab) and exit_status_of(tab) == 0, "no golden: table, status 0")
    tab, _ = cmp(case_a, now({"s1": 10.0, "s2": 10.0}), golden({"s1": 10.0}, name="other"))
    check(all(r[6] == "NEW" for r in tab), "case not in the golden: NEW")
    tab, _ = cmp(case_a, now({"s1": 10.0, "s2": 10.0}), g + golden({"s1": 1.0}, name="gone"), names={"a"})
    check(any(r[0] == "gone" and r[6] == "WARN" for r in tab), "a golden case cases.toml no longer has: WARN")
    check("ratio" in format_table(tab) and "summary" in summary_line(tab, 1, []), "table and summary render")
    # thresholds from spread: a noisy golden tolerates +50 %, a tight one does not
    noisy = golden({"s1": 10.0, "s2": 10.0}, repeats=5, spread=0.1)   # (max-min)/median = 0.2 -> fail 1.8
    tab, _ = cmp(case_a, now({"s1": 15.0, "s2": 10.0}), noisy)
    check(exit_status_of(tab) == 0 and any(r[1] == "s1" and r[6] == "WARN" and r[5].endswith("*") for r in tab),
          f"spread-derived threshold: +50 % is a WARN, not a FAIL: {tab}")
    tab, _ = cmp(case_a, now({"s1": 19.0, "s2": 10.0}), noisy)
    check(exit_status_of(tab) == 1, "...but +90 % fails")
    tab, _ = cmp(case_a, now({"s1": 4.0, "s2": 10.0}), g)
    check(_ == {"FAST"}, "FAST sets the re-record hint")

    # golden round trip: formats, order, partial merge
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "g", "x.csv")
        rows = golden({"s2": 1.5, "s1": 0.25}, e=-75.123456789012, it=11)
        rows[0]["gradient_rms"] = 0.0123
        order = {"z": 0, "a": 1}
        write_golden(path, rows + golden({"q": 2.0}, name="z"), order)
        back = read_golden(path)
        check([r["case"] for r in back] == ["z", "a", "a"] and [r["stage"] for r in back] == ["q", "s1", "s2"],
              f"rows in cases.toml order, then stage: {[(r['case'], r['stage']) for r in back]}")
        with open(path) as fh:
            lines = fh.read().splitlines()
        check(lines[0] == ",".join(GOLDEN_FIELDS), "golden header")
        check(",1.500000," in lines[3] and "-75.1234567890" in lines[3] and "1.2300000000e-02" in lines[3], f"formats: {lines[3]}")
        check(back[2]["median_s"] == 1.5 and back[2]["iterations"] == 11 and back[2]["repeats"] == 3
              and back[2]["gradient_rms"] == 0.0123 and back[1]["gradient_rms"] is None and back[2]["arch"] == "p",
              f"golden round trip {back[2]}")
        check(read_golden(os.path.join(tmp, "absent.csv")) == [], "absent golden reads as empty")
        newr = golden({"s1": 3.0, "s3": 1.0})
        ch = describe_changes(rows, newr)
        check(any("changed" in x and "s1" in x for x in ch) and any("new" in x and "s3" in x for x in ch)
              and any("removed" in x and "s2" in x for x in ch), f"change description {ch}")
        other = [dict(r, case="z", host="otherhost") for r in rows]
        merged = merge_golden(other + rows, newr)
        check({r["case"] for r in merged} == {"a", "z"} and len([r for r in merged if r["case"] == "a"]) == 2,
              "merge keeps other cases, replaces the ones that ran")
        check(all(r["host"] == "otherhost" for r in merged if r["case"] == "z"), "kept rows keep their own provenance")

    # build facts and the mismatch check
    info = parse_build_line("terco version 1\nbuild: device=yes gpu_arch=cc70 lmax=3 df_lmax=2 commit=abc compiler=NVHPC 26.5\n")
    check(info == {"compiler": "NVHPC 26.5", "device": "yes", "gpu_arch": "cc70", "lmax": "3", "df_lmax": "2",
                   "commit": "abc"}, f"{info}")
    check(parse_build_line("terco version 1\n") is None, "no build line")
    b_dev = "device=yes;gpu_arch=cc70;lmax=3;df_lmax=2"
    check(build_mismatch(b_dev, "NVHPC 26.5", b_dev, "NVHPC 25.1") == [], "same build, other compiler version: fine")
    d = build_mismatch(b_dev, "NVHPC 26.5", "device=yes;gpu_arch=cc80;lmax=2;df_lmax=2", "GNU 15.1.0")
    check(len(d) == 3 and any("gpu_arch: golden cc70" in x for x in d) and any("lmax: golden 3, this binary 2" in x for x in d)
          and any("compiler family" in x for x in d), f"mismatch sentences {d}")
    check(any(x.startswith("device:") for x in build_mismatch("device=no;lmax=2", "GNU 15", "device=yes;lmax=2", "GNU 15")),
          "device differs")
    check(build_mismatch("device=no;lmax=2", "GNU 15", "device=no;lmax=2;df_lmax_aux=4", "GNU 15") == [],
          "a key the golden does not carry is not compared")
    check(build_string({"device": False, "device_word": "no", "lmax": "2", "gpu_arch": "cc70", "df_lmax": "2"})
          == "device=no;lmax=2;df_lmax=2", "a host build string has no gpu_arch")

    # GPU selection and slugs
    def refuses(*a):
        try:
            resolve_gpu(*a)
        except UsageError:
            return True
        return False

    check(resolve_gpu(False, None, False, {}) is None, "host build needs no GPU")
    check(refuses(True, None, False, {}), "GPU build with no GPU chosen refuses")
    check(refuses(True, 0, False, {}), "--gpu 0 refuses")
    check(refuses(True, -1, False, {}), "--gpu -1 refuses")
    check(resolve_gpu(True, 0, True, {}) == "0", "--gpu 0 with --allow-gpu0")
    check(resolve_gpu(True, 2, False, {}) == "2", "--gpu 2")
    check(resolve_gpu(True, None, False, {"CUDA_VISIBLE_DEVICES": "3"}) == "3", "env GPU accepted")
    check(resolve_gpu(True, None, False, {"CUDA_VISIBLE_DEVICES": "GPU-1234-abcd"}) == "GPU-1234-abcd"
          and gpu_is_uuid("GPU-1234-abcd") and not gpu_is_uuid("2"), "a GPU-<uuid> in the environment is allowed")
    check(refuses(True, None, False, {"CUDA_VISIBLE_DEVICES": "0"}), "env GPU 0 refuses")
    check(refuses(True, None, False, {"CUDA_VISIBLE_DEVICES": "0,1"}), "env starting at GPU 0 refuses")
    check(refuses(True, None, False, {"CUDA_VISIBLE_DEVICES": ""}), "empty env refuses")
    check(refuses(True, None, False, {"CUDA_VISIBLE_DEVICES": "-1"}), "an unrecognisable env refuses")
    for name, slug in (("Tesla V100-DGXS-32GB", "v100"), ("Tesla V100-SXM2-32GB", "v100"), ("NVIDIA A100-SXM4-80GB", "a100"),
                       ("NVIDIA H100 80GB HBM3", "h100"), ("NVIDIA H200", "h200"),
                       ("NVIDIA GH200 480GB", "gh200"), ("NVIDIA Grace Hopper Superchip", "gh200"),
                       ("NVIDIA L40S", "l40s"), ("NVIDIA GeForce RTX 4090", "rtx-4090")):
        check(gpu_slug(name) == slug, f"slug {name} -> {gpu_slug(name)} (want {slug})")
    check(cpu_slug("model name\t: Intel(R) Xeon(R) Gold 6248 CPU @ 2.50GHz\n") == "cpu-intel-xeon-gold-6248",
          cpu_slug("model name\t: Intel(R) Xeon(R) Gold 6248 CPU @ 2.50GHz\n"))
    check(cpu_model_name("processor\t: 0\nmodel name\t: Foo CPU\n") == "Foo CPU", "cpu model name")
    ws = machine_warnings([dict(r, gpu_name="Tesla V100-SXM2-32GB", cpu_model="A", host="h1") for r in golden({"s1": 1.0})],
                          [case_a], {"gpu_name": "Tesla V100-DGXS-32GB", "cpu_model": "A"})
    check(len(ws) == 1 and "V100-SXM2" in ws[0], f"same family, other board: warned once {ws}")

    # the shipped cases.toml parses and every case is well formed
    defaults_real, _, cases_real = load_config(os.path.join(HERE, "cases.toml"))
    check(len(cases_real) >= 10, "cases.toml has the cases")
    used = set()
    for cs in cases_real:
        check(os.path.isfile(cs["geometry_path"]), f"geometry of {cs['name']} exists")
        check(isinstance(cs["args"], list) and cs["stages"], f"{cs['name']} has args and stages")
        check(cs["judge"], f"{cs['name']} judges at least one stage")
        check("total_wall" not in cs["judge"] and "scf.wall" not in cs["judge"], f"{cs['name']} does not judge a noisy total")
        used.add(os.path.basename(cs["geometry"]))
    systems = {f for f in os.listdir(os.path.join(HERE, "systems")) if f.endswith(".xyz")}
    check(systems == used, f"every geometry in systems/ is used by a case: unused {systems - used}")
    names = " ".join(cs["name"] for cs in cases_real)
    check(any("6-31g" == cs["args"][1] for cs in cases_real if cs["args"][:1] == ["-b"]), "a 6-31G case")
    check(any("6-31g*" == cs["args"][1] for cs in cases_real if cs["args"][:1] == ["-b"]), "a 6-31G* case")
    check("ecp" in names, "an ECP case")
    check(any("--no-incremental" not in cs["args"] and cs["name"].startswith("fock") for cs in cases_real),
          "a case on the default incremental path")

    # end to end with a fake terco: record, compare, failures, errors
    real_capture = run_capture
    with tempfile.TemporaryDirectory() as tmp:
        def make_fake(name, device, lmax, compiler="GNU 15.1.0"):
            path = os.path.join(tmp, name)
            with open(path, "w") as fh:
                fh.write(FAKE_TERCO.format(python=sys.executable, device=device, lmax=lmax, compiler=compiler,
                                           gpu_arch="gpu_arch=cc70 " if device == "yes" else ""))
            os.chmod(path, 0o755)
            return path

        fake = make_fake("terco_host", "no", 2)
        fake_gpu = make_fake("terco_gpu", "yes", 3, "NVHPC 26.5")
        cfg = os.path.join(tmp, "cases.toml")
        with open(cfg, "w") as fh:
            fh.write(FAKE_CASES)
        open(os.path.join(tmp, "geom.xyz"), "w").write("1\n\nH 0 0 0\n")
        gold = os.path.join(tmp, "goldens")
        counter = os.path.join(tmp, "counter")

        def run(extra, env=None, binary=None, repeats="2", out="out", capture=False, cwd=None):
            e = dict(os.environ)
            e.update(env or {})
            argv = ["--terco", binary or fake, "--config", cfg, "--goldens", gold, "--repeats", repeats] + extra
            if out:
                argv += ["--out", out if os.path.isabs(out) else os.path.join(tmp, out)] if cwd is None else ["--out", out]
            old_out, old_err, old_cwd = sys.stdout, sys.stderr, os.getcwd()
            sink = tempfile.TemporaryFile("w+")
            sys.stdout = sink
            sys.stderr = open(os.devnull, "w")
            try:
                if cwd:
                    os.chdir(cwd)
                try:
                    code = main(argv, environ=e)
                except UsageError:
                    code = 2
                sink.seek(0)
                text = sink.read()
            finally:
                os.chdir(old_cwd)
                sys.stderr.close()
                sys.stderr = old_err
                sys.stdout = old_out
                sink.close()
            return (code, text) if capture else code

        check(run([]) == 0, "no golden: exit 0")
        check(not os.path.exists(gold) or not os.listdir(gold), "no golden is written without --record")
        check(run(["--record"], repeats="5") == 0, "record (5 repeats): exit 0")
        res_rows = read_golden(os.path.join(tmp, "out", "result.csv"))
        check(len(res_rows) > 0 and all(r["repeats"] == 5 for r in res_rows), "--out result.csv in the golden format")
        check(os.path.isfile(os.path.join(tmp, "out", "a.1.json")), "--out keeps the raw JSON of every repeat")
        files = os.listdir(gold)
        check(len(files) == 1 and files[0].startswith("cpu-"), f"host golden named cpu-*: {files}")
        gfile = os.path.join(gold, files[0])
        rows = read_golden(gfile)
        check({r["case"] for r in rows} == {"a", "c"}, "case b (needs lmax3) was skipped on an lmax=2 build")
        check([r["case"] for r in rows] == sorted((r["case"] for r in rows)), "recorded rows follow cases.toml order")
        check(all(r["binary_commit"] == "abc1234" and r["omp_threads"] != "" for r in rows), "binary commit and thread provenance")
        check(run([], repeats="5") == 0, "second run against the golden: all OK")
        code, text = run([], {"FAKE_SCALE": "1.5"}, capture=True)
        check(code == 1 and "FAIL" in text, "a slowed binary exits 1")
        check(run([], {"FAKE_SCALE": "0.5"}) == 0, "a faster binary does not fail")
        code, text = run([], {"FAKE_SCALE": "0.5"}, capture=True)
        check("hint: " in text and "re-record" in text, "FAST prints the re-record hint")
        code, text = run([], {"FAKE_SCALE": "1.5"}, capture=True)
        fails = [ln for ln in text.splitlines() if " FAIL " in ln and "(run)" not in ln]
        check(len(fails) == 3, f"judged stages only (a: fock mean + gradient, c: fock mean), not total_wall/scf.wall: {len(fails)}")
        check(run([], {"FAKE_ESHIFT": "1e-5"}) == 1, "a wrong energy exits 1")
        check(run([], {"FAKE_ITERS": "12"}) == 0, "an iteration change alone is only a WARN")
        # one repeat among several carries a wrong energy
        if os.path.exists(counter):
            os.remove(counter)
        check(run([], {"FAKE_COUNTER": counter, "FAKE_ESHIFT_AT": "2"}, repeats="3") == 1, "one repeat's wrong energy exits 1")
        if os.path.exists(counter):
            os.remove(counter)
        code, text = run([], {"FAKE_COUNTER": counter, "FAKE_ITERS_AT": "2"}, repeats="3", capture=True)
        check("scf iterations repeats" in text, "repeats with different iteration counts are reported")
        # warm-up runs are discarded
        if os.path.exists(counter):
            os.remove(counter)
        check(run([], {"FAKE_COUNTER": counter, "FAKE_SLOWFIRST": "1"}, repeats="2") == 1,
              "without warm-up a slow first run (median of 2) shows")
        if os.path.exists(counter):
            os.remove(counter)
        check(run(["--warmup", "1"], {"FAKE_COUNTER": counter, "FAKE_SLOWFIRST": "1"}, repeats="2") == 0,
              "--warmup 1 discards the slow first run")
        check(int(open(counter).read()) == 2 * (1 + 2), "warm-up + repeats were all run (2 cases x 3 runs)")
        check(os.path.isfile(os.path.join(tmp, "out", "a.warmup1.json")), "warm-up output is kept apart")
        # an errored case does not abort the others
        code, text = run([], {"FAKE_FAIL_ARG": "z"}, capture=True)
        check(code == 2, f"an errored case: exit 2, got {code}")
        check(any(ln.startswith("c ") and "ERROR" in ln and ".log" in ln for ln in text.splitlines()),
              "the ERROR row names the log")
        check(any(ln.startswith("a ") and "OK" in ln for ln in text.splitlines()), "the other case is still judged")
        code, text = run([], {"FAKE_NOCONV": "1"}, capture=True)
        check(code == 2 and "NOT CONVERGED" in text, "exit status 2 with JSON is NOT CONVERGED")
        # a timeout kills the whole process group
        pidfile = os.path.join(tmp, "childpid")
        code, text = run([], {"FAKE_HANG": "1", "FAKE_SPAWN": pidfile}, capture=True)
        check(code == 2 and "ERROR" in text, "a timeout is an ERROR row, exit 2")
        alive = True
        try:
            pid = int(open(pidfile).read())
            for _ in range(50):
                try:
                    os.kill(pid, 0)
                    # a zombie still answers signal 0; check /proc state
                    with open(f"/proc/{pid}/stat") as fh:
                        if fh.read().split(")")[-1].split()[0] == "Z":
                            alive = False
                            break
                except (ProcessLookupError, FileNotFoundError):
                    alive = False
                    break
                import time
                time.sleep(0.1)
        except (OSError, ValueError):
            alive = True
        check(not alive, "the timeout killed the child's own child (process group)")
        check(run(["--cases", "nope"]) == 2, "unknown case exits 2")
        check(run(["--gpu", "-1"]) == 2, "--gpu -1 exits 2")
        check(run(["--dry-run"]) == 0 and run(["--list"]) == 0, "--dry-run and --list")
        # a relative --out works (the child runs in the output directory)
        reldir = os.path.join(tmp, "work")
        os.makedirs(reldir)
        check(run([], out="rel_out", cwd=reldir) == 0, "a relative --out runs")
        check(os.path.isfile(os.path.join(reldir, "rel_out", "a.1.json")), "...and keeps its files under the caller's cwd")
        # doctored goldens
        doctored = [dict(r, median_s=r["median_s"] / 2, min_s=r["min_s"] / 2, max_s=r["max_s"] / 2)
                    if r["stage"] == "scf.fock_mean_after_first" else r for r in rows]
        write_golden(gfile, doctored)
        check(run([], repeats="5") == 1, "a doctored golden (time halved) is caught: exit 1")
        wrong_e = [dict(r, energy_hartree=r["energy_hartree"] + 1e-4) for r in rows]
        write_golden(gfile, wrong_e)
        code, text = run([], repeats="5", capture=True)
        check(code == 1 and any("energy [Eh]" in ln and "FAIL" in ln for ln in text.splitlines()),
              "a wrong-energy golden is caught: exit 1")
        # a golden from another build is refused, per case, and not recorded into
        other_build = [dict(r, build="device=no;lmax=3;df_lmax=2") for r in rows]
        write_golden(gfile, other_build)
        code, text = run([], repeats="2", capture=True)
        check(code == 2 and "BUILD MISMATCH" in text and "lmax: golden 3, this binary 2" in text,
              f"build mismatch is refused and says what differs ({code})")
        before = open(gfile).read()
        check(run(["--record"], repeats="2") == 2 and open(gfile).read() == before, "--record into a mismatched golden is refused")
        check(run(["--allow-build-mismatch"], repeats="2") == 0, "--allow-build-mismatch compares anyway")
        other_cc = [dict(r, compiler="NVHPC 26.5") for r in rows]
        write_golden(gfile, other_cc)
        check(run([], repeats="2") == 2, "a different compiler family is refused")
        write_golden(gfile, rows)
        # a partial record keeps the other case's rows and provenance
        keep = [dict(r, host="elsewhere") if r["case"] == "c" else r for r in rows]
        write_golden(gfile, keep)
        check(run(["--record", "--cases", "a"], repeats="5") == 0, "partial record")
        after = read_golden(gfile)
        check({r["case"] for r in after} == {"a", "c"} and all(r["host"] == "elsewhere" for r in after if r["case"] == "c")
              and all(r["host"] != "elsewhere" for r in after if r["case"] == "a"),
              "a partial record keeps the other cases' rows and provenance")
        # GPU build: refuses without a GPU, never GPU 0 by default; identity from a canned nvidia-smi
        env_nogpu = {"CUDA_VISIBLE_DEVICES": ""}
        check(run([], env_nogpu, fake_gpu) == 2, "GPU build with no GPU chosen: exit 2")
        check(run(["--gpu", "0"], env_nogpu, fake_gpu) == 2, "GPU build on --gpu 0: exit 2")
        check(run(["--dry-run", "--gpu", "2"], env_nogpu, fake_gpu) == 0, "GPU build dry run on --gpu 2")
        seen_env = []

        def canned(cmd, env=None, timeout=60):
            seen_env.append(env.get("CUDA_DEVICE_ORDER") if env else None)
            if "--query-gpu=name,driver_version" in cmd:
                return 0, "Tesla V100-SXM2-32GB, 550.54.15\n", ""
            return 0, "| NVIDIA-SMI 550.54.15   Driver Version: 550.54.15   CUDA Version: 12.4 |\n", ""

        globals()["run_capture"] = canned
        try:
            gold_gpu = os.path.join(tmp, "goldens_gpu")
            argv = ["--terco", fake_gpu, "--config", cfg, "--goldens", gold_gpu, "--repeats", "1", "--gpu", "2",
                    "--record", "--out", os.path.join(tmp, "outgpu"), "--warmup", "0"]
            sink = open(os.devnull, "w")
            old = sys.stdout
            sys.stdout = sink
            try:
                code = main(argv, environ=dict(os.environ))
            finally:
                sys.stdout = old
                sink.close()
            check(code == 0 and os.path.isfile(os.path.join(gold_gpu, "v100.csv")), "device record: file named by the family")
            grows = read_golden(os.path.join(gold_gpu, "v100.csv"))
            check(grows and grows[0]["gpu_name"] == "Tesla V100-SXM2-32GB" and grows[0]["driver"] == "550.54.15"
                  and grows[0]["cuda"] == "12.4", "full GPU name, driver and CUDA version recorded")
            check(set(seen_env) == {"PCI_BUS_ID"}, "nvidia-smi is queried with CUDA_DEVICE_ORDER=PCI_BUS_ID")
            check("gpu_arch=cc70" in grows[0]["build"], "a device build records gpu_arch")
            check({r["case"] for r in grows} == {"a", "b", "c"}, "an lmax=3 device build also runs case b")
        finally:
            globals()["run_capture"] = real_capture

    print(f"self-test: {t.count} checks, {len(t.failures)} failed")
    return 1 if t.failures else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except UsageError as exc:
        print(f"error: {exc}", file=sys.stderr)
        sys.exit(2)
