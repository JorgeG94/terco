#!/usr/bin/env python3
"""The ONE place pyscf runs for terco's test suite.

Jorge's rule, from the scf_rks_vs_pyscf/UKS red build (pyscf's own UKS did
not reach conv_tol 1e-12 in 200 cycles -- terco agreed with it to 3.1e-11 at
terco's density, 4.1e-9 SCF to SCF, and the build still went red): pyscf
never runs in CI. It runs here, offline, on a developer machine with pyscf
installed, and the numbers it produces are committed as plain text under
this directory. Every `check_*_vs_pyscf`/`scf_*` Fortran program reads its
own file back and compares -- no pyscf, no numpy, in CI, ever.

Run with the pyscf virtualenv's interpreter, e.g.:

    /home/jorge/dev/metalquicha/.venv/bin/python3.12 test/pyscf_refs/generate.py

See README.md in this directory for the full story (which pyscf version,
the convergence-refusal rule, how to regenerate one case). Short version:

  * every FIXED-density comparison (Fock, erf-Fock, DF, gradient) runs the
    corresponding `check_*` binary from a CPU build to get terco's fixed,
    deterministic (seeded-LCG or closed-form) geometry/basis/density --
    never terco's SCF output -- then evaluates pyscf on exactly that input
    and commits pyscf's numbers, already permuted and rescaled into THIS
    BUILD's own AO order/normalisation, so the Fortran reader can diff
    straight across with no permutation of its own;
  * every SCF comparison (rks, scf_robust, rimp2) converges pyscf tightly
    (conv_tol 1e-12, more cycles, and a second-order solver where the
    first-order one does not reach it) and REFUSES to write a case whose
    SCF did not converge -- `SystemExit`, not a silently-wrong number;
  * every committed file carries a HASH of its own inputs (geometry, basis
    JSON, the Fortran source that builds the fixed input, this script) and
    the pyscf version, so check_staleness.py can tell a developer to rerun
    this script before trusting a stale number.
"""
from __future__ import annotations

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pyscf
from pyscf import df, gto, mp, scf, dft
from pyscf.dft import libxc, numint

HERE = Path(__file__).resolve().parent              # test/pyscf_refs
TEST = HERE.parent                                  # test/
ROOT = TEST.parent                                  # repo root
sys.path.insert(0, str(ROOT / "scripts"))
from cart_perm import ao_permutation  # noqa: E402

PYSCF_VERSION = pyscf.__version__

BUILD_DIR = Path(os.environ.get("TERCO_PYSCF_REF_BUILD", str(ROOT / "build-cpu")))
#: A CPU build (gfortran is fine; this never touches the GPU) whose test/
#: binaries this script runs to get a fixed geometry/basis/density -- see
#: this file's own header. Override with TERCO_PYSCF_REF_BUILD=<dir> to
#: point at a different one (e.g. the LMAX=3 build, for the f/g cases).

_SYMBOL = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 9: 'F', 10: 'Ne',
           14: 'Si', 15: 'P', 16: 'S', 17: 'Cl'}


# =========================================================================
# Shared AO-order / AO-normalisation machinery.
#
# MOVED HERE VERBATIM from the *_ref.py scripts this generator replaces
# (rys_grad_ref.py's bse_json_basis/terco_to_pyscf_perm, scf_rks_ref.py's
# inline (atom, l, exponent) permutation) rather than re-derived: both
# scripts read the identical BSE JSON convention and hit the identical
# pyscf cart=True AO-ordering quirk, and a second hand-copy of either is a
# second place for the two to silently drift apart.
# =========================================================================

def bse_json_basis(path, elements):
    """BSE JSON's electron_shells -> pyscf's basis={...} dict. Two shapes
    share this field:
      - len(angular_momentum) > 1: a COMBINED shell (6-31G-family "SP").
        One row of `coefficients` per l, same order as `angular_momentum`.
      - len(angular_momentum) == 1: a GENERAL CONTRACTION at that one l
        (cc-pVDZ oxygen's s-block: 9 exponents, 3 rows, all l=0). Every row
        is that SAME l.
    """
    d = json.load(open(path))
    out = {}
    for sym, z in elements.items():
        el = d['elements'][str(z)]
        shells = []
        for sh in el['electron_shells']:
            ls = sh['angular_momentum']
            exps = [float(e) for e in sh['exponents']]
            rows = sh['coefficients']
            if len(ls) > 1:
                pairs = zip(ls, rows)
            else:
                pairs = ((ls[0], row) for row in rows)
            for l, row in pairs:
                co = [float(c) for c in row]
                shells.append([l] + [[e, c] for e, c in zip(exps, co)])
        out[sym] = shells
    return out


_NCART = {0: 1, 1: 3, 2: 6, 3: 10, 4: 15}


def terco_to_pyscf_perm(atoms, basis_arg):
    """terco keeps every basis function in the LITERAL order its shell list
    gives them. pyscf's cart=True `Mole.build()` stable-sorts a given
    atom's shells by ascending angular momentum instead. Returns `perm`
    with `perm[i] = j` meaning terco's AO slot i is pyscf's AO slot j."""
    perm = []
    base = 0
    for sym, _ in atoms:
        shells = basis_arg[sym]
        ls = [row[0] for row in shells]
        n_at = sum(_NCART[l] for l in ls)
        terco_off = []
        off = 0
        for l in ls:
            terco_off.append(off)
            off += _NCART[l]
        order = sorted(range(len(ls)), key=lambda k: ls[k])
        pyscf_off = [0] * len(ls)
        off = 0
        for k in order:
            pyscf_off[k] = off
            off += _NCART[ls[k]]
        for k, l in enumerate(ls):
            for c in range(_NCART[l]):
                perm.append(base + pyscf_off[k] + c)
        base += n_at
    return np.array(perm, dtype=np.int64)


def ao_norm_from_overlap(s_terco, s_pyscf):
    """n_i = sqrt(S_terco(i,i)/S_pyscf(i,i)): the per-AO rescale between
    this build's normalisation and pyscf's cart=True one (their d/f
    off-diagonal Cartesian components are not all unit-normalised the same
    way). D_pyscf = N D_terco N (contravariant); G_pyscf = G_terco/(n_i n_j)
    for an ordinary <i|op|j> matrix like S, J, K, G (covariant) -- see
    rys_fock_ref.py's history in git log for the bug that confusing the two
    directions produces."""
    n = np.sqrt(np.diag(s_terco) / np.diag(s_pyscf))
    pred = n[:, None] * s_pyscf * n[None, :]
    err = np.abs(s_terco - pred).max()
    return n, err


# =========================================================================
# Hashing / staleness metadata, and the plain-text writer.
# =========================================================================

def sha256_of(*paths: Path) -> str:
    h = hashlib.sha256()
    for p in paths:
        h.update(Path(p).read_bytes())
    return h.hexdigest()


def inputs_hash(*paths) -> str:
    """The staleness fingerprint: every file whose bytes determine this
    case's fixed input or expected output -- geometry, basis JSON, the
    Fortran source that builds the density/grid, and this script itself.
    check_staleness.py (stdlib only, no pyscf) recomputes this and fails
    loudly if it no longer matches the HASH line in the committed file."""
    return sha256_of(Path(__file__), *[Path(p) for p in paths])


class RefWriter:
    """Accumulates the tagged, whitespace-delimited records
    pyscf_ref_io.F90 reads back (see that file's own header for the exact
    grammar) and writes them as one plain-text file."""

    def __init__(self, relname: str):
        self.path = HERE / relname
        self.lines = []

    def comment(self, text):
        self.lines.append("# " + text)

    def header(self, case, input_paths):
        rels = [str(Path(p).resolve().relative_to(ROOT)) for p in
                [Path(__file__), *[Path(p) for p in input_paths]]]
        self.comment("terco pyscf reference v1")
        self.comment("case: " + case)
        self.comment("generator: test/pyscf_refs/generate.py")
        self.comment("pyscf: " + PYSCF_VERSION)
        self.comment("hash: " + inputs_hash(*input_paths))
        self.comment("files: " + " ".join(rels))
        #: check_staleness.py (stdlib only) recomputes the hash from exactly
        #: this file list -- it never hardcodes per-case knowledge of what
        #: feeds into which reference, so a case that starts reading a new
        #: input is automatically covered the next time this is run.

    def int_(self, tag, v):
        self.lines.append(f"{tag} {int(v)}")

    def real(self, tag, v):
        self.lines.append(f"{tag} {float(v):.17e}")

    def str_(self, tag, v):
        self.lines.append(f"{tag} {v}")

    def _chunked(self, flat, perline=6):
        for i in range(0, len(flat), perline):
            self.lines.append(" ".join(f"{x:.17e}" for x in flat[i:i + perline]))

    def vector(self, tag, v):
        v = np.asarray(v, dtype=np.float64).ravel(order="F")
        self.lines.append(f"{tag} {len(v)}")
        self._chunked(v)

    def matrix(self, tag, a):
        a = np.asarray(a, dtype=np.float64)
        assert a.ndim == 2
        self.lines.append(f"{tag} {a.shape[0]} {a.shape[1]}")
        self._chunked(a.flatten(order="F"))

    def matrix3(self, tag, a):
        a = np.asarray(a, dtype=np.float64)
        assert a.ndim == 3
        self.lines.append(f"{tag} {a.shape[0]} {a.shape[1]} {a.shape[2]}")
        self._chunked(a.flatten(order="F"))

    def end(self):
        self.lines.append("END")

    def write(self):
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text("\n".join(self.lines) + "\n")
        print(f"  wrote {self.path.relative_to(ROOT)}")


class ConvergenceError(RuntimeError):
    """pyscf's SCF did not converge for a case this generator was about to
    commit. Refuse to write it -- see this file's own header."""


def converge_tight(make_mf, dm0, label, conv_tol=1.0e-12, conv_tol_grad=1.0e-8, max_cycle=200):
    """Converge an SCF to BOTH an energy and a GRADIENT tolerance, refuse if
    it cannot, and verify the result with a fresh mf before trusting it.

    pyscf's own convergence test is |dE| < conv_tol AND |g| < conv_tol_grad
    (default conv_tol_grad = sqrt(conv_tol) ~ 1e-6 at conv_tol=1e-12) --
    loose enough that a slow case (OH/wB97X-V UKS, measured) can stop up to
    ~1e-8 above its true minimum while still reporting `converged`. That is
    exactly what test/scf_rks.F90's own at-density control exposed: terco's
    energy AT the cached density agreed with pyscf to 2.5e-11 on every
    compiler, but the CACHED energy itself was the under-converged one by
    ~1.2e-8, past the 1e-8 gate the Fortran side holds the comparison to.
    An SCF energy is variational -- under-convergence only leaves it too
    HIGH -- so a comparison this tight needs `conv_tol_grad` explicit and
    tight, not pyscf's default.

    `make_mf()` must return a FRESH mf object (new instance, same mol/xc/
    grids) each call: once to do the actual convergence, once after it, at
    the SAME converged density, to cross-check `energy_tot(dm)` reproduces
    `e_tot` from an mf that never iterated at all -- catching a stale
    intermediate (cached J/K, grids,...) on the FIRST mf that the energy
    alone would not. Falls back to `.newton()` (second-order) when the
    first-order DIIS solver does not reach `conv_tol_grad` -- OH/PBE-UKS and
    OH/wB97X-V-UKS have both needed it, measured.

    Returns (mf, grad_norm). Raises ConvergenceError (never returns an
    unconverged or unverified result) otherwise.
    """
    mf = make_mf()
    mf.conv_tol = conv_tol
    mf.conv_tol_grad = conv_tol_grad
    mf.max_cycle = max_cycle
    mf.verbose = 0
    mf.kernel(dm0=dm0)
    if not mf.converged:
        mf2 = mf.newton()
        mf2.conv_tol = conv_tol
        mf2.conv_tol_grad = conv_tol_grad
        mf2.max_cycle = max_cycle
        mf2.kernel(dm0=dm0)
        mf = mf2
    if not mf.converged:
        raise ConvergenceError(
            f"{label}: pyscf SCF did not reach conv_tol={conv_tol}, "
            f"conv_tol_grad={conv_tol_grad} (DIIS and .newton() both tried) "
            f"-- refusing to write a reference from an under-converged "
            f"energy.")
    dm_final = mf.make_rdm1()
    # get_grad's third positional argument is a FOCK matrix, not a density
    # -- omit it and let pyscf build both from mo_coeff/mo_occ itself,
    # consistent with dm_final (mf.make_rdm1() with no arguments uses the
    # same mo_coeff/mo_occ).
    grad = mf.get_grad(mf.mo_coeff, mf.mo_occ)
    grad_norm = float(np.linalg.norm(grad))
    if grad_norm > conv_tol_grad:
        # Belt and suspenders: pyscf's own `converged` flag already implies
        # this, but the achieved norm is recorded in the reference header
        # regardless, so a future tightening of conv_tol_grad has a number
        # to check against instead of only a pass/fail.
        raise ConvergenceError(
            f"{label}: converged=True but |g|={grad_norm:.3e} > "
            f"conv_tol_grad={conv_tol_grad:.3e} -- pyscf's own check and "
            f"this one disagree; investigate before trusting either.")
    mf_fresh = make_mf()
    mf_fresh.verbose = 0
    e_fresh = mf_fresh.energy_tot(dm_final)
    if abs(e_fresh - mf.e_tot) > 1.0e-11:
        raise ConvergenceError(
            f"{label}: self-check failed -- a FRESH mf's energy_tot(dm) at "
            f"the converged density ({e_fresh:.14f}) disagrees with the "
            f"converged mf's own e_tot ({mf.e_tot:.14f}) by "
            f"{abs(e_fresh - mf.e_tot):.3e}, past 1e-11. Something the "
            f"converged mf cached (grids, J/K, ...) is not reproduced by a "
            f"plain energy_tot call -- investigate before trusting e_tot.")
    return mf, grad_norm


def run_bin(relexe, args, env_extra=None):
    """Run a test/-built binary from BUILD_DIR/test, with cwd=test/ (so its
    relative file reads -- water.xyz, basis_sets/... -- resolve exactly as
    they do under ctest)."""
    exe = BUILD_DIR / "test" / relexe
    if not exe.exists():
        raise FileNotFoundError(
            f"{exe} does not exist -- build a CPU config first (see "
            f"README.md), or point TERCO_PYSCF_REF_BUILD at one that has it.")
    env = dict(os.environ)
    env["TERCO_CART_ORDER"] = "lexicographic"
    #: The references this script commits are for TERCO_CART_ORDER=
    #: lexicographic ONLY -- the CI default (CMakeLists.txt's own default,
    #: and CI never overrides it). A gamess-order build's within-shell
    #: checks (check_mult/check_ao/check_xc_energy's pyscf halves) are a
    #: known gap against such a build -- see README.md.
    if env_extra:
        env.update(env_extra)
    # Not check=True: a check binary that now ALSO compares against a
    # reference this script has not written yet (first run, or a case
    # whose reference is being regenerated) exits 1 after writing its
    # dump -- the dump is what this call is for, and the caller verifies
    # it was written.
    proc = subprocess.run([str(exe)] + [str(a) for a in args], cwd=str(TEST), env=env)
    if proc.returncode != 0:
        print(f"  (note: {relexe} exited {proc.returncode} -- fine here if it already wrote its dump; "
              f"this run only needs the dump, not the comparison)")


def read_stream(path, dtype, count):
    return np.fromfile(path, dtype=dtype, count=count)


class BinReader:
    """Sequential reader over a terco `access='stream'` dump, the same
    layout the old `*_ref.py` scripts read."""

    def __init__(self, path):
        self.f = open(path, "rb")

    def i8(self, n=1):
        v = np.fromfile(self.f, dtype=np.int64, count=n)
        return int(v[0]) if n == 1 else v

    def f8(self, n):
        return np.fromfile(self.f, dtype=np.float64, count=n)

    def mat(self, n, m=None):
        m = n if m is None else m
        return self.f8(n * m).reshape(n, m, order="F")

    def mat3(self, n1, n2, n3):
        return self.f8(n1 * n2 * n3).reshape(n1, n2, n3, order="F")

    def str_(self, n):
        return self.f.read(n).decode().strip()

    def close(self):
        self.f.close()


if TERCO_CART_ORDER := os.environ.get("TERCO_CART_ORDER"):
    if TERCO_CART_ORDER not in ("lexicographic", "libcint"):
        sys.exit(f"generate.py must run with TERCO_CART_ORDER unset or "
                  f"'lexicographic' (got {TERCO_CART_ORDER!r}) -- see this "
                  f"file's run_bin() for why.")


# =========================================================================
# check_mult / check_ao / check_xc_energy's shared uncontracted probe
# system: s/p/d, one primitive each, three off-axis centres. Defined once
# here, in the SAME literal numbers check_mult.F90/check_ao.F90/
# check_xc_energy.F90 carry (that is what makes it reproducible without
# running the binary for check_mult; check_ao and check_xc_energy still
# need their own grid/density dumps, below).
# =========================================================================

def _probe_mol():
    coords = [(0.0, 0.0, 0.0), (1.3, 0.0, 0.7), (-0.4, 1.1, 1.9)]
    labels = ['Be', 'Li', 'He']
    sh_l, bas, atoms = [], {}, []
    for i, (lab, c) in enumerate(zip(labels, coords), start=1):
        bas[lab] = [[0, [0.9 + 0.17 * i, 1.0]],
                    [1, [1.4 + 0.11 * i, 1.0]],
                    [2, [1.9 + 0.13 * i, 1.0]]]
        atoms.append([lab, c])
        sh_l += [0, 1, 2]
    mol = gto.M(atom=atoms, basis=bas, cart=True, unit='Bohr', spin=None, verbose=0)
    return mol, sh_l


def gen_mult():
    """check_mult.F90's 39 Cartesian multipole components (through the
    octupole) against pyscf's int1e_r/int1e_rr/int1e_rrr. No Fortran binary
    needed -- the probe system is fully fixed (see _probe_mol); this is a
    direct port of mult_ref.py."""
    mol, sh_l = _probe_mol()
    nao = mol.nao_nr()
    s = mol.intor('int1e_ovlp_cart')
    r1 = mol.intor('int1e_r_cart').reshape(3, nao, nao)
    r2 = mol.intor('int1e_rr_cart').reshape(9, nao, nao)
    r3 = mol.intor('int1e_rrr_cart').reshape(27, nao, nao)
    ref = np.concatenate([r1, r2, r3], axis=0)   # (39, nao, nao), pyscf-lex order

    perm = ao_permutation(sh_l)   # pyscf-lex -> this build's order
    s = s[np.ix_(perm, perm)]
    ref = ref[:, perm][:, :, perm]

    w = RefWriter("mult.ref")
    w.header("check_mult", [ROOT / "test" / "check_mult.F90"])
    w.int_("NAO", nao)
    w.int_("NMULT", ref.shape[0])
    w.matrix("S", s)
    w.matrix3("MM", np.moveaxis(ref, 0, 2))   # (nao, nao, 39), matching trc_multipoles' layout
    w.end()
    w.write()


def gen_ao():
    """check_ao.F90's AO values/gradients at a probe cloud against pyscf's
    eval_ao. Needs the actual probe points (an LCG sequence) and terco's
    own chi/gchi, so it reads check_ao's own ao_probe.bin dump rather than
    re-deriving the LCG in Python -- see this module's own header."""
    mol, sh_l = _probe_mol()
    nao = mol.nao_nr()
    perm = ao_permutation(sh_l)

    run_bin("check_ao", [])
    dump = TEST / "ao_probe.bin"
    r = BinReader(dump)
    npts, n = r.i8(), r.i8()
    pts = r.f8(3 * npts).reshape(3, npts, order='F').T
    chi = r.f8(npts * n).reshape(npts, n, order='F')
    gchi = r.f8(3 * npts * n).reshape(3, npts, n, order='F')
    r.close()
    if n != nao:
        raise RuntimeError(f"gen_ao: nao mismatch terco {n} pyscf {nao}")

    ao = numint.eval_ao(mol, np.ascontiguousarray(pts), deriv=1)   # (4, npts, nao), pyscf-lex
    ao = ao[..., perm]

    w = RefWriter("ao.ref")
    w.header("check_ao", [ROOT / "test" / "check_ao.F90"])
    # Not `dump`: the probe bin is ephemeral (deleted at the end of this
    # function) and a function of check_ao.F90/this script alone, both
    # already hashed -- see RefWriter.header's own comment.
    w.int_("NPTS", npts)
    w.int_("NAO", nao)
    w.matrix("CHI", ao[0])
    w.matrix3("GCHI", np.moveaxis(ao[1:4], 0, 2))   # (npts, nao, 3)
    w.end()
    w.write()
    dump.unlink()


def gen_xc_energy():
    """check_xc_energy.F90's E_xc/V_xc, RKS and UKS, against pyscf's numint
    on the SAME grid and density (the seeded-LCG C C^T construction) --
    reads xc_probe.bin for the grid/density/terco-shape, same reasoning as
    gen_ao above."""
    mol, sh_l = _probe_mol()
    nao = mol.nao_nr()
    perm = ao_permutation(sh_l)
    inv = np.argsort(perm)

    NAMES = ['lda_x', 'svwn', 'pbe', 'blyp', 'b3lyp', 'b3lyp5', 'wb97x-v', 'wb97x', 'wb97',
             'vwn5', 'pw91', 'bp86', 'rpbe', 'o3lyp', 'revb3lyp', 'x3lyp', 'bhandh', 'sogga11x',
             'bhhlyp', 'hflyp']
    PYSCF = {'lda_x': 'LDA_X', 'svwn': 'LDA_X,LDA_C_VWN', 'pbe': 'PBE', 'blyp': 'BLYP',
             'b3lyp': 'HYB_GGA_XC_B3LYP', 'b3lyp5': 'HYB_GGA_XC_B3LYP5',
             'wb97x-v': 'HYB_GGA_XC_WB97X_V', 'wb97x': 'HYB_GGA_XC_WB97X', 'wb97': 'HYB_GGA_XC_WB97',
             'vwn5': 'LDA_C_VWN_RPA', 'pw91': 'GGA_X_PW91,GGA_C_PW91', 'bp86': 'GGA_X_B88,GGA_C_P86',
             'rpbe': 'GGA_X_RPBE,GGA_C_PBE', 'o3lyp': 'O3LYP', 'revb3lyp': 'REVB3LYP', 'x3lyp': 'X3LYP',
             'bhandh': '0.5*HF + 0.5*LDA_X, GGA_C_LYP', 'sogga11x': 'HYB_GGA_X_SOGGA11_X,GGA_C_SOGGA11_X',
             'bhhlyp': 'HYB_GGA_XC_BHANDHLYP', 'hflyp': 'HF,GGA_C_LYP'}

    run_bin("check_xc_energy", [])
    dump = TEST / "xc_probe.bin"
    r = BinReader(dump)
    npts, n, nfunc = r.i8(), r.i8(), r.i8()
    pts = r.f8(3 * npts).reshape(3, npts, order='F').T
    wts = r.f8(npts)
    dm = r.mat(n)
    results = []
    for _ in range(nfunc):
        exc, nelec = r.f8(2)
        vxc = r.mat(n)
        results.append((exc, nelec, vxc))
    dm2 = r.f8(2 * n * n).reshape(n, n, 2, order='F')
    uks = []
    for _ in range(nfunc):
        exc, nelec = r.f8(2)
        vxc = r.mat(n)
        uks.append((exc, nelec, vxc))
    r.close()
    if n != nao:
        raise RuntimeError(f"gen_xc_energy: nao mismatch terco {n} pyscf {nao}")

    grids = dft.gen_grid.Grids(mol)
    grids.coords = np.ascontiguousarray(pts)
    grids.weights = wts
    grids.non0tab = None
    ni = numint.NumInt()

    dm_lex = dm[np.ix_(inv, inv)]
    dm2_lex = (dm2[:, :, 0][np.ix_(inv, inv)], dm2[:, :, 1][np.ix_(inv, inv)])

    w = RefWriter("xc_energy.ref")
    w.header("check_xc_energy", [ROOT / "test" / "check_xc_energy.F90"])
    w.int_("NAO", nao)
    w.int_("NPTS", npts)
    w.int_("NFUNC", nfunc)
    for name in NAMES[:nfunc]:
        nelec_p, exc_p, vxc_p = ni.nr_rks(mol, grids, PYSCF[name], dm_lex)
        vxc_p = vxc_p[np.ix_(perm, perm)]
        w.real(f"EXC_{name}", exc_p)
        w.real(f"N_{name}", nelec_p)
        w.matrix(f"VXC_{name}", vxc_p)
        print(f"  xc_energy RKS {name:8s}  E_xc {exc_p:18.12f}  N {nelec_p:14.10f}")
    for name in NAMES[:nfunc]:
        nelec_p, exc_p, vxc_p = ni.nr_uks(mol, grids, PYSCF[name], dm2_lex)
        vxc_p = (vxc_p[0][np.ix_(perm, perm)], vxc_p[1][np.ix_(perm, perm)])
        w.real(f"UEXC_{name}", exc_p)
        w.real(f"UN_{name}", nelec_p.sum())
        w.matrix3(f"UVXC_{name}", np.stack(vxc_p, axis=2))
        print(f"  xc_energy UKS {name:8s}  E_xc {exc_p:18.12f}  N {nelec_p.sum():14.10f}")
    w.end()
    w.write()
    dump.unlink()


# =========================================================================
# VV10, the fixed-density, fixed-grid kernel check.
# =========================================================================

def gen_vv10():
    """check_vv10_vs_pyscf.F90's 3 cases (water/6-31G* RKS, water/cc-pVDZ
    RKS, OH/cc-pVDZ UKS) against pyscf's numint.nr_nlc_vxc, same grid,
    weights and density. See git log for the nlc_coeff provenance and the
    ADD-contract reasoning (vmat is seeded nonzero; what must match pyscf
    is vmat_after - vmat_seed)."""
    XC_NLC = 'WB97X_V'
    b_ref, c_ref = libxc.nlc_coeff(XC_NLC)[0][0]
    assert abs(b_ref - 6.0) < 1e-12 and abs(c_ref - 0.01) < 1e-12, \
        "wB97X-V NLC parameters drifted from (6.0, 0.01)"

    run_bin("check_vv10_vs_pyscf", [])
    dump = TEST / "vv10_probe.bin"
    bases = ['basis_sets/6-31gs.json', 'basis_sets/cc-pvdz.json', 'basis_sets/cc-pvdz.json']
    r = BinReader(dump)
    ncase = r.i8()
    w = RefWriter("vv10.ref")
    w.header("check_vv10_vs_pyscf",
             [ROOT / "test" / "check_vv10_vs_pyscf.F90"] +
             [TEST / b for b in sorted(set(bases))])
    w.int_("NCASE", ncase)
    for ic in range(ncase):
        kind = r.i8()
        nat, nao, npts = r.i8(), r.i8(), r.i8()
        at_z = r.i8(nat)
        rcoord = r.f8(3 * nat).reshape(3, nat, order='F').T
        coords = r.f8(3 * npts).reshape(3, npts, order='F').T
        wts = r.f8(npts)
        s_terco = r.mat(nao)
        if kind == 0:
            dmat_terco = r.mat(nao)
        else:
            da_terco = r.mat(nao)
            db_terco = r.mat(nao)
            dmat_terco = da_terco + db_terco
        if kind == 0:
            r.mat(nao)   # vmat0, unused (only the ADD contract on the Fortran side needs it)
        else:
            r.mat(nao); r.mat(nao)
        r.f8(1)          # e_nlc_terco, unused here
        if kind == 0:
            r.mat(nao)   # vmat_terco, unused
        else:
            r.mat(nao); r.mat(nao)

        elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
        basis_arg = bse_json_basis(TEST / bases[ic], elements)
        atoms = [(_SYMBOL[int(z)], tuple(rcoord[i])) for i, z in enumerate(at_z)]
        spin = 1 if kind == 1 else 0
        mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True, verbose=0, spin=spin)
        if mol.nao_nr() != nao:
            raise RuntimeError(f"gen_vv10 case {ic}: nao mismatch terco {nao} pyscf {mol.nao_nr()}")

        perm = terco_to_pyscf_perm(atoms, basis_arg)
        inv_perm = np.argsort(perm)
        s_terco_p = s_terco[np.ix_(inv_perm, inv_perm)]
        dmat_terco_p = dmat_terco[np.ix_(inv_perm, inv_perm)]
        s_pyscf = mol.intor('int1e_ovlp')
        n, s_err = ao_norm_from_overlap(s_terco_p, s_pyscf)
        if s_err >= 1e-9:
            raise RuntimeError(f"gen_vv10 case {ic}: overlap check failed, {s_err:.3e}")
        dmat_pyscf = n[:, None] * dmat_terco_p * n[None, :]

        grids = dft.gen_grid.Grids(mol)
        grids.coords = np.ascontiguousarray(coords)
        grids.weights = wts
        grids.non0tab = None
        ni = numint.NumInt()
        nelec_p, exc_p, vmat_p = ni.nr_nlc_vxc(mol, grids, XC_NLC, dmat_pyscf)
        # V_nlc is covariant like S (dE/dD): predicted-terco = N vmat_p N,
        # in terco's AO order.
        v_pred_terco = (n[:, None] * vmat_p * n[None, :])[np.ix_(perm, perm)]

        w.int_(f"KIND_{ic}", kind)
        w.real(f"E_NLC_{ic}", exc_p)
        w.matrix(f"VMAT_{ic}", v_pred_terco)
        print(f"  vv10 case {ic} ({'RKS' if kind == 0 else 'UKS'}): E_nlc(pyscf) = {exc_p:18.14f}")
    w.end()
    w.write()
    r.close()
    dump.unlink()


# =========================================================================
# Fixed-density four-centre checks: rys_fock, erf_fock, df_vs_pyscf.
# =========================================================================

def gen_rys_fock(case, xyz, json_rel, tol):
    dump = TEST / f"rys_fock_probe_{case}.bin"
    run_bin("check_rys_fock_vs_pyscf_" + case, [xyz, json_rel, dump.name])
    if not dump.exists():
        print(f"  rys_fock {case}: binary self-skipped (LMAX too low for this build); not regenerating")
        return
    r = BinReader(dump)
    nat, nao = r.i8(), r.i8()
    at_z = r.i8(nat)
    rcoord = r.f8(3 * nat).reshape(3, nat, order='F').T
    d_terco = r.mat(nao)
    g_terco = r.mat(nao)
    s_terco = r.mat(nao)
    r.close()

    elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
    basis_arg = bse_json_basis(TEST / json_rel, elements)
    atoms = [(_SYMBOL[int(z)], tuple(rcoord[i])) for i, z in enumerate(at_z)]
    mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True, verbose=0)
    if mol.nao_nr() != nao:
        raise RuntimeError(f"gen_rys_fock {case}: nao mismatch terco {nao} pyscf {mol.nao_nr()}")

    perm = terco_to_pyscf_perm(atoms, basis_arg)
    inv_perm = np.argsort(perm)
    s_terco_p = s_terco[np.ix_(inv_perm, inv_perm)]
    d_terco_p = d_terco[np.ix_(inv_perm, inv_perm)]

    s_pyscf = mol.intor('int1e_ovlp')
    n, s_err = ao_norm_from_overlap(s_terco_p, s_pyscf)
    if s_err >= tol:
        raise RuntimeError(f"gen_rys_fock {case}: overlap check failed, {s_err:.3e}")
    d_pyscf = n[:, None] * d_terco_p * n[None, :]

    eri = mol.intor('int2e', aosym='s1').reshape(nao, nao, nao, nao)
    j = np.einsum('ijkl,lk->ij', eri, d_pyscf)
    k = np.einsum('ijkl,jk->il', eri, d_pyscf)
    g_ref_pyscf = j - 0.5 * k
    g_ref_terco = (g_ref_pyscf * n[:, None] * n[None, :])[np.ix_(perm, perm)]

    w = RefWriter(f"rys_fock_{case}.ref")
    w.header("check_rys_fock_vs_pyscf_" + case,
             [ROOT / "test" / "check_rys_fock_vs_pyscf.F90", TEST / json_rel, TEST / xyz])
    w.int_("NAT", nat)
    w.int_("NAO", nao)
    w.matrix("S", s_terco)
    w.matrix("G", g_ref_terco)
    w.end()
    w.write()
    worst = np.abs(g_terco - g_ref_terco).max()
    print(f"  rys_fock {case}: max|G_terco - G_pyscf| (at generation time) = {worst:.3e}")
    dump.unlink()


def gen_erf_fock(case, xyz, json_rel):
    """omega=0.3 range-separated J/K against pyscf's own
    mol.with_range_coulomb -- one reference per (geometry, basis), reused
    by all three Fortran kernel paths (hgp/nosym/rys) since the fixed
    density and the expected pyscf J/K do not depend on which terco code
    path is under test, only on the basis and omega."""
    OMEGA = 0.3
    dump = TEST / f"erf_fock_probe_{case}.bin"
    # Any of the three kernel-path binaries (hgp/nosym/rys) builds the
    # identical fixed density and basis -- "_hgp" is simply the one that
    # always exists (the "_rys" variant only when TERCO_RYS is on).
    run_bin(f"check_erf_fock_vs_pyscf_{case}_hgp", [xyz, json_rel, dump.name, "0"])
    if not dump.exists():
        print(f"  erf_fock {case}: binary self-skipped (LMAX too low for this build); not regenerating")
        return
    r = BinReader(dump)
    nat, nao = r.i8(), r.i8()
    at_z = r.i8(nat)
    rcoord = r.f8(3 * nat).reshape(3, nat, order='F').T
    d_terco = r.mat(nao)
    r.mat(nao); r.mat(nao)   # gj_terco, gk_terco: not needed, every path recomputes its own
    s_terco = r.mat(nao)
    r.close()

    elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
    basis_arg = bse_json_basis(TEST / json_rel, elements)
    atoms = [(_SYMBOL[int(z)], tuple(rcoord[i])) for i, z in enumerate(at_z)]
    mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True, verbose=0)
    if mol.nao_nr() != nao:
        raise RuntimeError(f"gen_erf_fock {case}: nao mismatch terco {nao} pyscf {mol.nao_nr()}")

    perm = terco_to_pyscf_perm(atoms, basis_arg)
    inv_perm = np.argsort(perm)
    s_terco_p = s_terco[np.ix_(inv_perm, inv_perm)]
    d_terco_p = d_terco[np.ix_(inv_perm, inv_perm)]
    s_pyscf = mol.intor('int1e_ovlp')
    n, s_err = ao_norm_from_overlap(s_terco_p, s_pyscf)
    if s_err >= 1e-9:
        raise RuntimeError(f"gen_erf_fock {case}: overlap check failed, {s_err:.3e}")
    d_pyscf = n[:, None] * d_terco_p * n[None, :]

    with mol.with_range_coulomb(OMEGA):
        vj_ref, vk_ref = scf.hf.get_jk(mol, d_pyscf)
    j_terco = (vj_ref * n[:, None] * n[None, :])[np.ix_(perm, perm)]
    k_terco = (vk_ref * n[:, None] * n[None, :])[np.ix_(perm, perm)]

    w = RefWriter(f"erf_fock_{case}.ref")
    w.header("check_erf_fock_vs_pyscf_" + case,
             [ROOT / "test" / "check_erf_fock_vs_pyscf.F90", TEST / json_rel, TEST / xyz])
    w.int_("NAT", nat)
    w.int_("NAO", nao)
    w.matrix("S", s_terco)
    w.matrix("J", j_terco)
    w.matrix("K", k_terco)
    w.end()
    w.write()
    dump.unlink()


def gen_df_vs_pyscf():
    """trc_df_2c/trc_df_3c (water/cc-pVTZ orbital, cc-pVTZ-RI auxiliary --
    the first comparison with a g shell) against pyscf's
    int2c2e/aux_e2(int3c2e)."""
    xyz, json_rel, aux_rel = "water.xyz", "basis_sets/cc-pvtz.json", "basis_sets/cc-pvtz-rifit.json"
    dump = TEST / "df_vs_pyscf_probe.bin"
    run_bin("check_df_vs_pyscf", [xyz, json_rel, aux_rel, dump.name])
    if not dump.exists():
        print("  df_vs_pyscf: binary self-skipped (DF_LMAX/DF_LMAX_AUX too low for this build); not regenerating")
        return
    r = BinReader(dump)
    nat, nao, naux = r.i8(), r.i8(), r.i8()
    at_z = r.i8(nat)
    rcoord = r.f8(3 * nat).reshape(3, nat, order='F').T
    o_terco = r.mat(nao)
    j_terco = r.mat(naux)
    tens_terco = r.mat3(nao, nao, naux)
    r.close()

    elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
    basis_arg = bse_json_basis(TEST / json_rel, elements)
    aux_arg = bse_json_basis(TEST / aux_rel, elements)
    atoms = [(_SYMBOL[int(z)], tuple(rcoord[i])) for i, z in enumerate(at_z)]
    mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True, verbose=0)
    auxmol = gto.M(atom=atoms, basis=aux_arg, unit='Bohr', cart=True, verbose=0)
    if mol.nao_nr() != nao or auxmol.nao_nr() != naux:
        raise RuntimeError("gen_df_vs_pyscf: basis mismatch")

    perm_orb = terco_to_pyscf_perm(atoms, basis_arg)
    perm_aux = terco_to_pyscf_perm(atoms, aux_arg)
    inv_orb = np.argsort(perm_orb)
    inv_aux = np.argsort(perm_aux)
    o_terco_p = o_terco[np.ix_(inv_orb, inv_orb)]
    j_terco_p = j_terco[np.ix_(inv_aux, inv_aux)]

    o_pyscf = mol.intor('int2c2e')
    n_orb, o_err = ao_norm_from_overlap(o_terco_p, o_pyscf)
    if o_err >= 1e-10:
        raise RuntimeError(f"gen_df_vs_pyscf: orbital (i|j) metric check failed, {o_err:.3e}")
    j_pyscf = auxmol.intor('int2c2e')
    n_aux, j_err = ao_norm_from_overlap(j_terco_p, j_pyscf)
    if j_err >= 1e-10:
        raise RuntimeError(f"gen_df_vs_pyscf: 2c (P|Q) check failed, {j_err:.3e}")

    tens_pyscf = df.incore.aux_e2(mol, auxmol, 'int3c2e')
    tens_ref_terco = tens_pyscf * n_orb[:, None, None] * n_orb[None, :, None] * n_aux[None, None, :]
    tens_ref_terco = tens_ref_terco[np.ix_(perm_orb, perm_orb, perm_aux)]
    j_ref_terco = (j_pyscf * n_aux[:, None] * n_aux[None, :])[np.ix_(perm_aux, perm_aux)]

    w = RefWriter("df_vs_pyscf.ref")
    w.header("check_df_vs_pyscf", [ROOT / "test" / "check_df_vs_pyscf.F90",
                                    TEST / json_rel, TEST / aux_rel, TEST / xyz])
    w.int_("NAT", nat)
    w.int_("NAO", nao)
    w.int_("NAUX", naux)
    w.matrix("O", o_terco)
    w.matrix("J", j_ref_terco)
    w.matrix3("TENS", tens_ref_terco)
    w.end()
    w.write()
    worst = np.abs(tens_terco - tens_ref_terco).max()
    print(f"  df_vs_pyscf: max|tens_terco - tens_pyscf| (at generation time) = {worst:.3e}")
    dump.unlink()


def gen_rys_grad(case, json_rel):
    """eri%gradient against pyscf's int2e_ip1, contracted by hand (pyscf's
    own Gradients class assumes a converged SCF density; this is a fixed,
    non-SCF one, same reasoning as check_rys_grad_fd.F90). RHF and UHF, on
    water."""
    xyz = "water.xyz"
    dump = TEST / f"rys_grad_probe_{case}.bin"
    run_bin(f"check_rys_grad_vs_pyscf_{case}", [json_rel, dump.name])
    r = BinReader(dump)
    nat, nao = r.i8(), r.i8()
    at_z = r.i8(nat)
    rcoord = r.f8(3 * nat).reshape(3, nat, order='F').T
    drhf = r.mat(nao)
    da = r.mat(nao)
    db = r.mat(nao)
    grad_rhf_terco = r.f8(3 * nat).reshape(3, nat, order='F').T
    grad_uhf_terco = r.f8(3 * nat).reshape(3, nat, order='F').T
    s_terco = r.mat(nao)
    r.close()

    elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
    basis_arg = bse_json_basis(TEST / json_rel, elements)
    atoms = [(_SYMBOL[int(z)], tuple(rcoord[i])) for i, z in enumerate(at_z)]
    mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True, verbose=0)
    if mol.nao_nr() != nao:
        raise RuntimeError(f"gen_rys_grad {case}: nao mismatch terco {nao} pyscf {mol.nao_nr()}")

    perm = terco_to_pyscf_perm(atoms, basis_arg)
    inv_perm = np.argsort(perm)
    s_terco_p = s_terco[np.ix_(inv_perm, inv_perm)]
    drhf_p = drhf[np.ix_(inv_perm, inv_perm)]
    da_p = da[np.ix_(inv_perm, inv_perm)]
    db_p = db[np.ix_(inv_perm, inv_perm)]

    s_pyscf = mol.intor('int1e_ovlp')
    n, s_err = ao_norm_from_overlap(s_terco_p, s_pyscf)
    if s_err >= 1e-9:
        raise RuntimeError(f"gen_rys_grad {case}: overlap check failed, {s_err:.3e}")
    drhf_pyscf = n[:, None] * drhf_p * n[None, :]
    da_pyscf = n[:, None] * da_p * n[None, :]
    db_pyscf = n[:, None] * db_p * n[None, :]

    eri1 = mol.intor('int2e_ip1', comp=3).reshape(3, nao, nao, nao, nao)
    aoslice = mol.aoslice_by_atom()

    def grad2e_rhf(dm):
        vj = -np.einsum('xijkl,lk->xij', eri1, dm)
        vk = -np.einsum('xijkl,jk->xil', eri1, dm)
        vhf = vj - 0.5 * vk
        de = np.zeros((nat, 3))
        for ia in range(nat):
            p0, p1 = aoslice[ia, 2:]
            de[ia] = np.einsum('xij,ij->x', vhf[:, p0:p1], dm[p0:p1]) * 4.0
        return de

    def grad2e_uhf(da_, db_):
        dtot = da_ + db_
        vj_tot = -np.einsum('xijkl,lk->xij', eri1, dtot)
        vk_a = -np.einsum('xijkl,jk->xil', eri1, da_)
        vk_b = -np.einsum('xijkl,jk->xil', eri1, db_)
        veff_a = vj_tot - vk_a
        veff_b = vj_tot - vk_b
        de = np.zeros((nat, 3))
        for ia in range(nat):
            p0, p1 = aoslice[ia, 2:]
            de[ia] = (np.einsum('xij,ij->x', veff_a[:, p0:p1], da_[p0:p1]) +
                      np.einsum('xij,ij->x', veff_b[:, p0:p1], db_[p0:p1])) * 4.0
        return de

    grad_rhf_py = grad2e_rhf(drhf_pyscf)
    grad_uhf_py = grad2e_uhf(da_pyscf, db_pyscf)

    w = RefWriter(f"rys_grad_{case}.ref")
    w.header("check_rys_grad_vs_pyscf_" + case,
              [ROOT / "test" / "check_rys_grad_vs_pyscf.F90", TEST / json_rel, TEST / xyz])
    w.int_("NAT", nat)
    w.int_("NAO", nao)
    w.matrix("S", s_terco)
    w.matrix("GRAD_RHF", grad_rhf_py.T)   # (3, nat), matching terco's own layout
    w.matrix("GRAD_UHF", grad_uhf_py.T)
    w.end()
    w.write()
    worst_rhf = np.abs(grad_rhf_terco - grad_rhf_py).max()
    worst_uhf = np.abs(grad_uhf_terco - grad_uhf_py).max()
    print(f"  rys_grad {case}: RHF worst|terco-pyscf| {worst_rhf:.3e}  UHF {worst_uhf:.3e} (at generation time)")
    dump.unlink()


# =========================================================================
# SCF comparisons: converge pyscf tightly, REFUSE an unconverged case.
# =========================================================================

def gen_rks():
    """scf_rks.F90's 12 cases (water RHF/SVWN/PBE/BLYP/B3LYP/wB97X-V/
    wB97X/wB97 restricted; OH radical UHF/PBE/wB97X-V/wB97X unrestricted),
    6-31G uncontracted to primitives. Commits pyscf's own converged energy
    AND its own converged density (in terco's AO order) -- the at-density
    control is then: terco evaluates ITS energy at that committed density
    (trc_scf_run with max_iter=1 and dguess=this density; see
    src/trc_scf.F90's own comment '# --- energy at this density ---',
    which is exactly pyscf's energy_tot(dm0) equivalent -- a single
    iteration's etot is computed from the density the iteration STARTED
    with, before the Roothaan step updates it) and compares to the SAME
    committed energy."""
    SYMBOL = {1: 'H', 6: 'C', 7: 'N', 8: 'O'}
    PYSCF = {'svwn': 'LDA_X,LDA_C_VWN', 'pbe': 'PBE', 'blyp': 'BLYP', 'b3lyp': 'HYB_GGA_XC_B3LYP',
             'b3lyp5': 'HYB_GGA_XC_B3LYP5', 'lda_x': 'LDA_X', 'pbe0': 'PBE0',
             'wb97x-v': 'wb97x-v', 'wb97x': 'wb97x', 'wb97': 'wb97'}

    run_bin("scf_rks", [])
    dump = TEST / "rks_probe.bin"
    r = BinReader(dump)
    ncase = r.i8()
    w = RefWriter("rks.ref")
    w.header("scf_rks", [ROOT / "test" / "scf_rks.F90",
                          TEST / "basis_sets" / "6-31g.json", TEST / "water.xyz", TEST / "oh.xyz"])
    w.int_("NCASE", ncase)
    for ic in range(ncase):
        n, natm, nsh, nspin, na, nb = r.i8(), r.i8(), r.i8(), r.i8(), r.i8(), r.i8()
        z = r.f8(natm)
        rr = r.f8(3 * natm).reshape(3, natm, order='F').T
        shells = []
        for _ in range(nsh):
            a, l = r.i8(), r.i8()
            e = r.f8(1)[0]
            shells.append((int(a), int(l), float(e)))
        name = r.str_(8)
        npts = r.i8()
        if npts:
            pts = r.f8(3 * npts).reshape(3, npts, order='F').T
            wts = r.f8(npts)
        npts_nlc = r.i8()
        if npts_nlc:
            pts_nlc = r.f8(3 * npts_nlc).reshape(3, npts_nlc, order='F').T
            wts_nlc = r.f8(npts_nlc)
        e_t = r.f8(1)[0]
        dm = r.f8(n * n * nspin).reshape(n, n, nspin, order='F')

        atoms, bas = [], {}
        for a in range(natm):
            lab = f'{SYMBOL[int(round(z[a]))]}{a + 1}'
            atoms.append([lab, tuple(rr[a])])
            bas[lab] = [[l, [e, 1.0]] for (aa, l, e) in shells if aa == a + 1]
        charge = int(round(z.sum())) - (na + nb)
        mol = gto.M(atom=atoms, basis=bas, cart=True, unit='Bohr', charge=charge, spin=na - nb, verbose=0)
        if mol.nao_nr() != n:
            raise RuntimeError(f"gen_rks case {ic}: nao mismatch terco {n} pyscf {mol.nao_nr()}")

        off, k = {}, 0
        for (aa, l, e) in shells:
            off[(aa, l, round(e, 10))] = k
            k += (l + 1) * (l + 2) // 2
        perm = []
        for ib in range(mol.nbas):
            key = (mol.bas_atom(ib) + 1, mol.bas_angular(ib), round(float(mol.bas_exp(ib)[0]), 10))
            perm.extend(range(off[key], off[key] + mol.bas_len_cart(ib)))
        perm = np.array(perm)
        inv_perm = np.argsort(perm)
        dm_pyscf_order = dm[perm][:, perm]

        def make_mf():
            if name:
                m = dft.RKS(mol) if nspin == 1 else dft.UKS(mol)
                m.xc = PYSCF[name]
                m.grids.coords = np.ascontiguousarray(pts)
                m.grids.weights = wts
                m.grids.non0tab = None
                if npts_nlc:
                    m.nlcgrids.coords = np.ascontiguousarray(pts_nlc)
                    m.nlcgrids.weights = wts_nlc
                    m.nlcgrids.non0tab = None
            else:
                m = scf.RHF(mol) if nspin == 1 else scf.UHF(mol)
            return m

        label = (name or 'hf') + ('/U' if nspin == 2 else '/R')
        dm0 = dm_pyscf_order[:, :, 0] if nspin == 1 else (dm_pyscf_order[:, :, 0], dm_pyscf_order[:, :, 1])
        mf, grad_norm = converge_tight(make_mf, dm0, f"rks case {ic} ({label})")
        e_p = mf.e_tot

        dm_conv = mf.make_rdm1()
        if nspin == 1:
            dm_conv_terco = dm_conv[np.ix_(inv_perm, inv_perm)][:, :, None]
        else:
            dm_conv_terco = np.stack([dm_conv[0][np.ix_(inv_perm, inv_perm)],
                                      dm_conv[1][np.ix_(inv_perm, inv_perm)]], axis=2)

        w.str_(f"NAME_{ic}", name or "hf")
        w.int_(f"NSPIN_{ic}", nspin)
        w.int_(f"NAO_{ic}", n)
        w.real(f"E_{ic}", e_p)
        w.real(f"CONV_TOL_{ic}", mf.conv_tol)
        w.real(f"CONV_TOL_GRAD_{ic}", mf.conv_tol_grad)
        w.real(f"GRAD_NORM_{ic}", grad_norm)
        w.matrix3(f"DM_{ic}", dm_conv_terco)
        de = abs(e_t - e_p)
        print(f"  rks case {ic:2d} {label:9s} E(terco,old-live-compare) {e_t:18.12f}  "
              f"E(pyscf,converged) {e_p:18.12f}  |dE| {de:.2e}  |g| {grad_norm:.2e}")
    w.end()
    w.write()
    dump.unlink()


def gen_scf_robust():
    """scf_robust.F90's 4 hard cases: cache pyscf's own converged energy,
    at pyscf's own basis library where it has the name (6-31G, aug-cc-pVDZ)
    or the same BSE JSON terco reads (d-aug-cc-pVDZ, not in pyscf's
    library)."""
    SYMBOL = {1: 'H', 7: 'N', 8: 'O', 10: 'Ne'}
    run_bin("scf_robust", [])
    dump = TEST / "scf_robust_probe.bin"
    r = BinReader(dump)
    ncase = r.i8()
    w = RefWriter("scf_robust.ref")
    w.header("scf_robust", [ROOT / "test" / "scf_robust.F90"])
    w.int_("NCASE", ncase)
    for ic in range(ncase):
        n, natm, nspin, na, nb, charge = (r.i8() for _ in range(6))
        xyzname = r.str_(16)
        basname = r.str_(16)
        e_terco = r.f8(1)[0]

        atoms = []
        with open(TEST / (xyzname + ".xyz")) as g:
            nat = int(g.readline())
            g.readline()
            for _ in range(nat):
                parts = g.readline().split()
                atoms.append((parts[0], (float(parts[1]), float(parts[2]), float(parts[3]))))

        try:
            mol = gto.M(atom=atoms, basis=basname, unit='Angstrom', cart=True, charge=charge,
                        spin=na - nb, verbose=0)
        except Exception:
            elements = {sym: {'H': 1, 'N': 7, 'O': 8, 'Ne': 10}[sym] for sym, _ in atoms}
            bpath = TEST / "basis_sets" / (basname + ".json")
            basis_arg = bse_json_basis(bpath, elements)
            mol = gto.M(atom=atoms, basis=basis_arg, unit='Angstrom', cart=True, charge=charge,
                        spin=na - nb, verbose=0)
        if mol.nao_nr() != n:
            raise RuntimeError(f"gen_scf_robust {xyzname}/{basname}: nao mismatch terco {n} pyscf {mol.nao_nr()}")

        def make_mf():
            return scf.UHF(mol) if nspin == 2 else scf.RHF(mol)

        mf, grad_norm = converge_tight(make_mf, None, f"scf_robust {xyzname}/{basname}")
        e_py = mf.e_tot

        w.str_(f"XYZ_{ic}", xyzname)
        w.str_(f"BAS_{ic}", basname)
        w.real(f"E_{ic}", e_py)
        w.real(f"CONV_TOL_{ic}", mf.conv_tol)
        w.real(f"CONV_TOL_GRAD_{ic}", mf.conv_tol_grad)
        w.real(f"GRAD_NORM_{ic}", grad_norm)
        print(f"  scf_robust {xyzname:14s} {basname:14s}  E(terco,old) {e_terco:20.10f}  "
              f"E(pyscf,converged) {e_py:20.10f}  |diff| {abs(e_terco - e_py):10.3e}  |g| {grad_norm:.2e}")
    w.end()
    w.write()
    dump.unlink()


def _gen_rimp2_one(probe_relname, out_relname):
    SYMBOL = {1: 'H', 6: 'C', 7: 'N', 8: 'O'}
    dump = TEST / probe_relname
    r = BinReader(dump)
    n, natm, nsh, nsha, nocc = (r.i8() for _ in range(5))
    z = r.f8(natm)
    rr = r.f8(3 * natm).reshape(3, natm, order='F').T

    def shells(count):
        out = []
        for _ in range(count):
            a, l, npr = r.i8(), r.i8(), r.i8()
            e = r.f8(npr)
            c = r.f8(npr)
            out.append((a, l, [[float(x), float(y)] for x, y in zip(e, c)]))
        return out

    ao = shells(nsh)
    ax = shells(nsha)
    e_scf, e_os, e_ss = r.f8(3)
    r.close()

    atoms, bas, auxbas = [], {}, {}
    for a in range(natm):
        lab = f'{SYMBOL[int(round(z[a]))]}{a + 1}'
        atoms.append([lab, tuple(rr[a])])
        bas[lab] = [[l] + prims for (aa, l, prims) in ao if aa == a + 1]
        auxbas[lab] = [[l] + prims for (aa, l, prims) in ax if aa == a + 1]
    mol = gto.M(atom=atoms, basis=bas, cart=True, unit='Bohr', verbose=0)
    if mol.nao_nr() != n:
        raise RuntimeError(f"gen_rimp2 {probe_relname}: nao mismatch terco {n} pyscf {mol.nao_nr()}")
    def make_mf():
        return scf.RHF(mol)

    mf, grad_norm = converge_tight(make_mf, None, f"rimp2 {probe_relname} RHF")
    pt = mp.dfmp2.DFMP2(mf)
    pt.with_df = df.DF(mol, auxbasis=auxbas)
    pt.verbose = 0
    pt.kernel()
    ref_os, ref_ss = pt.e_corr_os, pt.e_corr_ss

    w = RefWriter(out_relname)
    w.header("scf_rimp2:" + probe_relname, [ROOT / "test" / "scf_rimp2.F90"])
    w.real("E_SCF", mf.e_tot)
    w.real("E_OS", ref_os)
    w.real("E_SS", ref_ss)
    w.real("CONV_TOL", mf.conv_tol)
    w.real("CONV_TOL_GRAD", mf.conv_tol_grad)
    w.real("GRAD_NORM", grad_norm)
    w.end()
    w.write()
    print(f"  rimp2 {probe_relname}: E_RHF(terco,old) {e_scf:20.12f}  E_RHF(pyscf) {mf.e_tot:20.12f}  "
          f"|diff| {abs(e_scf - mf.e_tot):9.2e}  |g| {grad_norm:.2e}")
    print(f"             E_os(terco,old) {e_os:20.12f}  E_os(pyscf) {ref_os:20.12f}  "
          f"|diff| {abs(e_os - ref_os):9.2e}")
    print(f"             E_ss(terco,old) {e_ss:20.12f}  E_ss(pyscf) {ref_ss:20.12f}  "
          f"|diff| {abs(e_ss - ref_ss):9.2e}")
    dump.unlink()


def gen_rimp2():
    """RI-MP2 on water (two cases: terco's synthetic primitive 6-31G/aux
    ladder, and the MolSSI BSE cc-pVDZ/cc-pVDZ-RI pair) against pyscf's
    RHF + DF-MP2 on the identical orbital and auxiliary sets."""
    run_bin("scf_rimp2", [])
    _gen_rimp2_one("rimp2_probe.bin", "rimp2.ref")
    run_bin("scf_rimp2", ["water.xyz", "basis_sets/cc-pvdz.json", "basis_sets/cc-pvdz-rifit.json",
                          "rimp2_ccpvdz_probe.bin"])
    _gen_rimp2_one("rimp2_ccpvdz_probe.bin", "rimp2_ccpvdz.ref")


# =========================================================================
# The translated XC kernels against libxc (folded in from the old
# test/xc_ref.py, unchanged logic, writing test/reference_xc.dat exactly
# as before -- check_xc.F90 needs no change, and never imported pyscf
# itself).
# =========================================================================

def gen_xc_kernels():
    sys.path.insert(0, str(ROOT / "tools"))
    from xc_kernel_discovery import DISAGREES_WITH_LIBXC  # noqa: E402

    KERNELS = {
        "slater_exchange": ("LDA_X", "lda"), "vwn": ("LDA_C_VWN", "lda"),
        "vwn_rpa": ("LDA_C_VWN_RPA", "lda"), "pbe_x": ("GGA_X_PBE", "gga"),
        "pbe_c": ("GGA_C_PBE", "gga"), "b88": ("GGA_X_B88", "gga"), "lyp": ("GGA_C_LYP", "gga"),
        "wb97x_v_xc": ("WB97X_V", "gga"), "wb97x_xc": ("WB97X", "gga"), "wb97_xc": ("WB97", "gga"),
        "wb97x_d_xc": ("WB97X_D", "gga"), "wb97x_d3_xc": ("WB97X_D3", "gga"),
        "vwn3": ("LDA_C_VWN_3", "lda"), "pz81": ("LDA_C_PZ", "lda"), "pz81_mod": ("LDA_C_PZ_MOD", "lda"),
        "pw91_lda": ("LDA_C_PW", "lda"), "pw91_lda_mod": ("LDA_C_PW_MOD", "lda"),
        "pw91_lda_rpa": ("LDA_C_PW_RPA", "lda"), "rev_pbe_x": ("GGA_X_PBE_R", "gga"),
        "rpbe_x": ("GGA_X_RPBE", "gga"), "optx_x": ("GGA_X_OPTX", "gga"),
        "mpw91_x": ("GGA_X_MPW91", "gga"), "pw86_x": ("GGA_X_PW86", "gga"), "pw91_x": ("GGA_X_PW91", "gga"),
        "p86_c": ("GGA_C_P86", "gga"), "p86vwn_ft_c": ("GGA_C_P86VWN_FT", "gga"),
        "pw91_c": ("GGA_C_PW91", "gga"), "pbe_sol_c": ("GGA_C_PBE_SOL", "gga"),
        "bmk_c": ("GGA_C_BMK", "gga"), "n12_c": ("GGA_C_N12", "gga"), "n12_sx_c": ("GGA_C_N12_SX", "gga"),
        "sogga11_x_c": ("GGA_C_SOGGA11_X", "gga"), "sogga11_x_x": ("HYB_GGA_X_SOGGA11_X", "gga"),
        "b97_d": ("GGA_XC_B97_D", "gga"),
        "ityh_x": ("GGA_X_ITYH", "gga", 0.2), "ityh_x_015": ("GGA_X_ITYH", "gga", 0.15),
        "ityh_x_033": ("GGA_X_ITYH", "gga", 0.33),
        "hjs_pbe_x": ("GGA_X_HJS_PBE", "gga", 0.11), "lrcwpbe_hjs_pbe_x": ("GGA_X_HJS_PBE", "gga", 0.3),
        "lrcwpbeh_hjs_pbe_x": ("GGA_X_HJS_PBE", "gga", 0.2),
    }

    RHO = np.logspace(-9, 1, 11)
    S_RED = np.array([0.0, 1e-3, 0.1, 1.0, 10.0, 1e3])
    ZETA = np.array([0.0, 0.3, 0.7, 0.99])
    COS_AB = np.array([-0.8, 0.0, 0.6])
    _SWITCH_DELTA = np.array([-1e-3, -1e-6, 0.0, 1e-6, 1e-3])

    def _rho_switch(omega):
        c = (9.0 ** (2.0 / 3.0)) * (0.682784063255295503581 ** 2) * omega \
            * 1.442249570307408382321638310780109588392 / 18.0
        return (c / 1.35) ** 3

    def wb97_extra_points(omega, polar):
        rho_switch = _rho_switch(omega)
        if not polar:
            pts = [(rho, s * rho ** (8.0 / 3.0))
                   for rho in rho_switch * (1.0 + _SWITCH_DELTA) for s in S_RED]
            return np.array(pts)
        rho_a_switch = rho_switch / 2.0
        pts = []
        other = 1.0e-3
        for ra_is_switch in (True, False):
            for d in _SWITCH_DELTA:
                ra = rho_a_switch * (1.0 + d) if ra_is_switch else other
                rb = other if ra_is_switch else rho_a_switch * (1.0 + d)
                for s in S_RED[::2]:
                    saa = s * ra ** (8.0 / 3.0)
                    sbb = s * rb ** (8.0 / 3.0)
                    for c in COS_AB:
                        pts.append((ra, rb, saa, c * np.sqrt(saa * sbb), sbb))
        return np.array(pts)

    _OMEGA = {"wb97x_v_xc": 0.3, "wb97x_xc": 0.3, "wb97_xc": 0.4,
              "wb97x_d_xc": 0.2, "wb97x_d3_xc": 0.25}

    def unpolar_points(family):
        pts = []
        for rho in RHO:
            if family == "lda":
                pts.append((rho,))
            else:
                for s in S_RED:
                    pts.append((rho, s * rho ** (8.0 / 3.0)))
        return np.array(pts)

    def polar_points(family):
        pts = []
        for rho in RHO:
            for z in ZETA:
                ra, rb = 0.5 * rho * (1 + z), 0.5 * rho * (1 - z)
                if family == "lda":
                    pts.append((ra, rb))
                else:
                    for s in S_RED[::2]:
                        saa = s * ra ** (8.0 / 3.0)
                        sbb = s * rb ** (8.0 / 3.0)
                        for c in COS_AB:
                            pts.append((ra, rb, saa, c * np.sqrt(saa * sbb), sbb))
        return np.array(pts)

    def libxc_rho(family, polar, pts):
        n = len(pts)
        if not polar:
            if family == "lda":
                return pts[:, 0]
            r = np.zeros((4, n))
            r[0] = pts[:, 0]
            r[1] = np.sqrt(pts[:, 1])
            return r
        if family == "lda":
            return np.stack([pts[:, 0], pts[:, 1]])
        ra, rb, saa, sab, sbb = pts.T
        a = np.zeros((4, n))
        b = np.zeros((4, n))
        a[0], b[0] = ra, rb
        a[1] = np.sqrt(saa)
        gb = np.sqrt(sbb)
        with np.errstate(invalid="ignore", divide="ignore"):
            cos = np.where(a[1] * gb > 0, sab / (a[1] * gb), 0.0)
        cos = np.clip(cos, -1.0, 1.0)
        b[1] = gb * cos
        b[2] = gb * np.sqrt(1.0 - cos ** 2)
        return np.stack([a, b])

    out = TEST / "reference_xc.dat"
    xc_inputs = [Path(__file__), ROOT / "tools" / "xc_kernel_discovery.py"]
    rels = [str(p.resolve().relative_to(ROOT)) for p in xc_inputs]
    lines = ["# kernel family polar npts nin nout", "# libxc %s via pyscf" % libxc.libxc_version(),
             "# hash: " + sha256_of(*xc_inputs), "# files: " + " ".join(rels)]
    for stem, spec in KERNELS.items():
        if stem in DISAGREES_WITH_LIBXC:
            continue
        name, family = spec[0], spec[1]
        eval_omega = spec[2] if len(spec) > 2 else None
        for polar in (0, 1):
            pts = polar_points(family) if polar else unpolar_points(family)
            if stem in _OMEGA:
                pts = np.concatenate([pts, wb97_extra_points(_OMEGA[stem], polar == 1)])
            rho = libxc_rho(family, polar, pts)
            exc, vxc, fxc, _ = libxc.eval_xc(name, rho, spin=polar, deriv=2, omega=eval_omega)
            cols = [exc]
            nv = 2 if polar else 1
            ns = 3 if polar else 1
            cols.append(vxc[0].reshape(len(pts), nv))
            if family == "gga":
                cols.append(vxc[1].reshape(len(pts), ns))
            cols.append(fxc[0].reshape(len(pts), 3 if polar else 1))
            if family == "gga":
                cols.append(fxc[1].reshape(len(pts), 6 if polar else 1))
                cols.append(fxc[2].reshape(len(pts), 6 if polar else 1))
            ref = np.column_stack([c.reshape(len(pts), -1) for c in cols])
            lines.append(f"{stem} {family} {polar} {len(pts)} {pts.shape[1]} {ref.shape[1]}")
            for p, r_ in zip(pts, ref):
                lines.append(" ".join(f"{x:.17e}" for x in np.concatenate([p, r_])))
            print(f"  xc kernel {stem:16s} polar={polar} {len(pts):4d} points, {ref.shape[1]:2d} outputs")
    out.write_text("\n".join(lines) + "\n")
    print(f"  wrote {out.relative_to(ROOT)}")


CASES = {
    "mult": gen_mult,
    "ao": gen_ao,
    "xc_energy": gen_xc_energy,
    "vv10": gen_vv10,
    "rys_grad_631gs": lambda: gen_rys_grad("631gs", "basis_sets/6-31gs.json"),
    "rys_grad_ccpvdz": lambda: gen_rys_grad("ccpvdz", "basis_sets/cc-pvdz.json"),
    "rys_fock_water_ccpvdz": lambda: gen_rys_fock("water_ccpvdz", "water.xyz", "basis_sets/cc-pvdz.json", 1e-9),
    "rys_fock_water_ccpvtz": lambda: gen_rys_fock("water_ccpvtz", "water.xyz", "basis_sets/cc-pvtz.json", 1e-10),
    "rys_fock_formamide_def2tzvp": lambda: gen_rys_fock(
        "formamide_def2tzvp", "formamide.xyz", "basis_sets/def2-tzvp.json", 1e-9),
    "erf_fock_water_631gs": lambda: gen_erf_fock("water_631gs", "water.xyz", "basis_sets/6-31gs.json"),
    "erf_fock_water_ccpvdz": lambda: gen_erf_fock("water_ccpvdz", "water.xyz", "basis_sets/cc-pvdz.json"),
    "erf_fock_formamide_ccpvdz": lambda: gen_erf_fock("formamide_ccpvdz", "formamide.xyz", "basis_sets/cc-pvdz.json"),
    "df_vs_pyscf": gen_df_vs_pyscf,
    "rks": gen_rks,
    "scf_robust": gen_scf_robust,
    "rimp2": gen_rimp2,
    "xc_kernels": gen_xc_kernels,
}


def main():
    import argparse
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("cases", nargs="*", choices=list(CASES),
                    help="case names to (re)generate; default: all")
    args = p.parse_args()
    todo = args.cases or list(CASES)
    print(f"pyscf {PYSCF_VERSION}, build dir {BUILD_DIR}")
    failures = []
    for name in todo:
        print(f"\n== {name} ==")
        try:
            CASES[name]()
        except ConvergenceError as exc:
            print(f"  REFUSED: {exc}")
            failures.append(name)
        except FileNotFoundError as exc:
            print(f"  SKIPPED (binary not built in {BUILD_DIR}): {exc}")
    if failures:
        sys.exit(f"\n{len(failures)} case(s) refused (pyscf did not converge): {failures}")
    print("\nAll requested references written.")


if __name__ == "__main__":
    main()
