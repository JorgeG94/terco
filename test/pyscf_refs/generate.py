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
    JSON, the Fortran source that builds the fixed input) and of the code
    in this script that wrote it -- its CASES entry and every definition
    that entry reaches, not the whole file (refhash.py) -- plus the pyscf
    version, so check_staleness.py can tell a developer to rerun exactly
    the cases whose own code or inputs changed.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy.linalg
import pyscf
from pyscf import df, gto, mp, scf, dft
from pyscf.dft import libxc, numint
from pyscf.grad import rhf as pyscf_grad_rhf

HERE = Path(__file__).resolve().parent              # test/pyscf_refs
TEST = HERE.parent                                  # test/
ROOT = TEST.parent                                  # repo root
sys.path.insert(0, str(ROOT / "scripts"))
from cart_perm import ao_permutation  # noqa: E402
sys.path.insert(0, str(HERE))
import refhash  # noqa: E402

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


_NSPH = {0: 1, 1: 3, 2: 5, 3: 7, 4: 9}


def terco_to_pyscf_perm_sph(atoms, basis_arg):
    """Same shell-level reordering as `terco_to_pyscf_perm`, for a
    SPHERICAL basis (block sizes 2l+1 rather than ncart(l)). No per-AO
    normalisation rescale is needed on top of this one (unlike the
    Cartesian case's `ao_norm_from_overlap`): scripts/gen_sph.py's own
    self-test checks terco's cart2sph blocks (SPH_T) against
    `pyscf.gto.mole.cart2sph(l, normalized='sp')` -- the exact table pyscf
    itself uses -- to 1e-13, so a spherical AO here and pyscf's cart=False
    one of the same (atom, l, exponent) are already identically
    normalised."""
    perm = []
    base = 0
    for sym, _ in atoms:
        shells = basis_arg[sym]
        ls = [row[0] for row in shells]
        n_at = sum(_NSPH[l] for l in ls)
        terco_off = []
        off = 0
        for l in ls:
            terco_off.append(off)
            off += _NSPH[l]
        order = sorted(range(len(ls)), key=lambda k: ls[k])
        pyscf_off = [0] * len(ls)
        off = 0
        for k in order:
            pyscf_off[k] = off
            off += _NSPH[ls[k]]
        for k, l in enumerate(ls):
            for c in range(_NSPH[l]):
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

#: The CASES key main() is running; every header records it, and its hash
#: covers that entry's code (refhash.py), not this whole file.
_CURRENT_CASE = None


def ref_header_lines(input_paths):
    """The generator/hash/files lines of a reference's header, for the case
    main() is running. The hash (refhash.case_hash) covers the CASES entry,
    every definition in this file it reaches, and the input files --
    geometry, basis JSON, the Fortran source that builds the density/grid.
    check_staleness.py (stdlib only, no pyscf) recomputes it with the same
    refhash module and fails loudly if it no longer matches."""
    if _CURRENT_CASE is None:
        raise RuntimeError("ref_header_lines: no current case -- run through main()/CASES")
    gen = refhash.Generator.from_path(Path(__file__))
    paths = [Path(p).resolve() for p in input_paths]
    rels = [str(p.relative_to(ROOT)) for p in paths]
    return [f'generator: test/pyscf_refs/generate.py CASES["{_CURRENT_CASE}"] -> '
            f'{gen.function(_CURRENT_CASE)}',
            "hash-scheme: " + refhash.SCHEME,
            "hash: " + refhash.case_hash(gen, _CURRENT_CASE, paths),
            "files: " + " ".join(rels)]


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
        generator, *hashing = ref_header_lines(input_paths)
        self.comment("terco pyscf reference v1")
        self.comment("case: " + case)
        self.comment(generator)
        self.comment("pyscf: " + PYSCF_VERSION)
        for line in hashing:
            self.comment(line)
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


def require_converged(mf, label):
    """Refuse an mf whose own `converged` flag is unset. _pcm_converge's
    check between its passes; converge_tight does the rest of the file's."""
    if not mf.converged:
        raise ConvergenceError(
            f"{label}: pyscf SCF did not converge (conv_tol={mf.conv_tol}, "
            f"conv_tol_grad={mf.conv_tol_grad}, max_cycle={mf.max_cycle}) -- "
            f"refusing to write a reference from an unconverged energy.")


def converge_tight(make_mf, dm0, label, conv_tol=1.0e-12, conv_tol_grad=1.0e-8, max_cycle=200,
                    self_check_tol=1.0e-11):
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
    if abs(e_fresh - mf.e_tot) > self_check_tol:
        raise ConvergenceError(
            f"{label}: self-check failed -- a FRESH mf's energy_tot(dm) at "
            f"the converged density ({e_fresh:.14f}) disagrees with the "
            f"converged mf's own e_tot ({mf.e_tot:.14f}) by "
            f"{abs(e_fresh - mf.e_tot):.3e}, past {self_check_tol:.1e}. "
            f"Something the converged mf cached (grids, J/K, ...) is not "
            f"reproduced by a plain energy_tot call -- investigate before "
            f"trusting e_tot.")
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
    # identical fixed density and basis; "_hgp" is the one used here.
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


def gen_1e_grad(case, xyz, json_rel):
    """trc_1e_gradient's three terms (dS, dT, dV) against pyscf, each in
    isolation, at a fixed non-SCF density (check_1e_grad_vs_pyscf.F90's own
    header has the full reasoning for why three separate pyscf quantities,
    not one).

    THE "ip" SIGN. Every pyscf `int1e_ip*` integral is <nabla_r mu|op|nu>,
    the gradient of the BRA FUNCTION wrt the electron coordinate -- for a
    Gaussian centred at A, that is the NEGATIVE of d/dA (shifting the
    coordinate one way is shifting the centre the other). So
    d<mu|op|nu>/dA_mu = -int1e_ip<op>, confirmed independently (not taken on
    a pyscf docstring's word) by a central finite difference of
    Tr(W, S(R)) and Tr(D, T(R))/Tr(D, V(R)) at a displaced geometry, fixed
    W/D: agreed to 8 figures at h=1e-5.

    PER-ATOM ASSEMBLY. Standard pyscf technique (grad.rhf.grad_elec's own
    overlap/kinetic terms): slice the bra-derivative tensor to the AO rows
    belonging to one atom (mol.aoslice_by_atom()), add the transpose for
    the symmetric ket contribution (an atom's shells appear on BOTH sides
    of a symmetric matrix like S/T/V), and contract with the fixed
    matrix -- no aoslice needed for a result already per-atom.

    THE OPERATOR TERM. dV cannot be split into basis/Hellmann-Feynman this
    way without re-deriving pyscf's own int1e_iprinv convention, so dV is
    built from pyscf's OWN validated combined kinetic+nuclear hcore
    derivative (grad.rhf.Gradients.hcore_generator, which already bundles
    the nuclear-attraction basis derivative and the Hellmann-Feynman term
    per atom correctly -- that is what hcore_generator is FOR) minus the
    kinetic piece built the same way as dT above: dV = hcore_deriv(atom) -
    dT_only(atom). Confirmed against the same finite difference as dT,
    same figures.
    """
    dump = TEST / f"onee_grad_probe_{case}.bin"
    run_bin(f"check_1e_grad_vs_pyscf_{case}", [xyz, json_rel, dump.name, case])
    r = BinReader(dump)
    nat, nao = r.i8(), r.i8()
    at_z = r.i8(nat)
    rcoord = r.f8(3 * nat).reshape(3, nat, order='F').T
    wmat = r.mat(nao)
    dmat = r.mat(nao)
    gs_terco = r.f8(3 * nat).reshape(3, nat, order='F').T
    gt_terco = r.f8(3 * nat).reshape(3, nat, order='F').T
    gv_terco = r.f8(3 * nat).reshape(3, nat, order='F').T
    s_terco = r.mat(nao)
    r.close()

    elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
    basis_arg = bse_json_basis(TEST / json_rel, elements)
    atoms = [(_SYMBOL[int(z)], tuple(rcoord[i])) for i, z in enumerate(at_z)]
    mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True, verbose=0)
    if mol.nao_nr() != nao:
        raise RuntimeError(f"gen_1e_grad {case}: nao mismatch terco {nao} pyscf {mol.nao_nr()}")

    perm = terco_to_pyscf_perm(atoms, basis_arg)
    inv_perm = np.argsort(perm)
    s_terco_p = s_terco[np.ix_(inv_perm, inv_perm)]
    w_terco_p = wmat[np.ix_(inv_perm, inv_perm)]
    d_terco_p = dmat[np.ix_(inv_perm, inv_perm)]

    s_pyscf = mol.intor('int1e_ovlp')
    n, s_err = ao_norm_from_overlap(s_terco_p, s_pyscf)
    if s_err >= 1e-9:
        raise RuntimeError(f"gen_1e_grad {case}: overlap check failed, {s_err:.3e}")
    w_pyscf = n[:, None] * w_terco_p * n[None, :]
    d_pyscf = n[:, None] * d_terco_p * n[None, :]

    aoslices = mol.aoslice_by_atom()
    s1 = -mol.intor('int1e_ipovlp', comp=3)
    t1 = -mol.intor('int1e_ipkin', comp=3)
    mf = scf.RHF(mol)
    hcore_deriv = pyscf_grad_rhf.Gradients(mf).hcore_generator(mol)

    def place(tensor, atm):
        p0, p1 = aoslices[atm, 2:]
        out = np.zeros((3, nao, nao))
        out[:, p0:p1] = tensor[:, p0:p1]
        return out + out.transpose(0, 2, 1)

    grad_s_py = np.zeros((nat, 3))
    grad_t_py = np.zeros((nat, 3))
    grad_v_py = np.zeros((nat, 3))
    for ia in range(nat):
        gs_ia = place(s1, ia)
        gt_ia = place(t1, ia)
        gv_ia = hcore_deriv(ia) - gt_ia
        grad_s_py[ia] = -np.einsum('xij,ij->x', gs_ia, w_pyscf)
        grad_t_py[ia] = np.einsum('xij,ij->x', gt_ia, d_pyscf)
        grad_v_py[ia] = np.einsum('xij,ij->x', gv_ia, d_pyscf)

    w = RefWriter(f"onee_grad_{case}.ref")
    w.header("check_1e_grad_vs_pyscf_" + case,
              [ROOT / "test" / "check_1e_grad_vs_pyscf.F90", TEST / json_rel, TEST / xyz])
    w.int_("NAT", nat)
    w.int_("NAO", nao)
    w.matrix("S", s_terco)
    w.matrix("GRAD_S", grad_s_py.T)
    w.matrix("GRAD_T", grad_t_py.T)
    w.matrix("GRAD_V", grad_v_py.T)
    w.end()
    w.write()
    worst_s = np.abs(gs_terco - grad_s_py).max()
    worst_t = np.abs(gt_terco - grad_t_py).max()
    worst_v = np.abs(gv_terco - grad_v_py).max()
    print(f"  onee_grad {case}: worst|terco-pyscf| dS {worst_s:.3e}  dT {worst_t:.3e}  "
          f"dV {worst_v:.3e} (at generation time)")
    dump.unlink()


# =========================================================================
# SCF comparisons: converge pyscf tightly, REFUSE an unconverged case.
# =========================================================================

def gen_full_grad():
    """check_full_grad_vs_pyscf.F90's three cases: the FULL analytic
    RHF/UHF gradient (trc_scf_grad%trc_scf_gradient -- one-electron,
    eri%gradient's two-electron part, nuclear repulsion) against pyscf's
    own `Gradients` class AT A CONVERGED density, on SPHERICAL cc-pVDZ
    bases -- water/RHF, formamide/RHF, OH/UHF. This is the one comparison
    that exercises trc_scf_gradient's own spherical fold (D and W folded
    from nao_sph to the Cartesian nao the kernels need, T D T^T); every
    other gradient reference in this file (onee_grad_*, rys_grad_*) feeds
    `trc_1e_gradient`/`eri%gradient` a synthetic density directly, at the
    Cartesian nao, and never touches that fold.

    FIXED INITIAL DENSITY, water/formamide. NOT NEEDED: a closed-shell
    RHF ground state at this geometry has no other stationary point
    nearby, so terco (its own default GWH guess) and pyscf (its own
    default "minao" guess) converging independently land at the SAME
    density, same reasoning rks.ref itself uses for ITS restricted cases
    (which also commit no DM0 -- see that function's own header).

    FIXED INITIAL DENSITY, OH. NEEDED, same reasoning and the same recipe
    as rks.ref's own OH cases: OH's unpaired electron sits in one of a
    rotationally degenerate pair of pi orbitals, a genuine zero-curvature
    direction, so terco and pyscf's own independent default guesses are
    not guaranteed to pick the SAME orientation. A core-Hamiltonian guess
    perturbed by a substantial field along y breaks the degeneracy
    strongly enough to pick a definite one; that same density (DM0,
    already converged, permuted into terco's own spherical AO order) is
    what test/check_full_grad_vs_pyscf.F90 seeds trc_scf_run's `dguess`
    with, so both sides converge toward the SAME rotation rather than two
    that merely share an energy.
    """
    w = RefWriter("full_grad_ccpvdz.ref")
    w.header("check_full_grad_vs_pyscf",
             [ROOT / "test" / "check_full_grad_vs_pyscf.F90", TEST / "basis_sets" / "cc-pvdz.json",
              TEST / "water.xyz", TEST / "formamide.xyz", TEST / "oh.xyz"])
    cases = [
        ("water_ccpvdz", "water", 0, 1),
        ("formamide_ccpvdz", "formamide", 0, 1),
        ("oh_ccpvdz", "oh", 0, 2),
    ]
    for name, xyzname, charge, mult in cases:
        atoms = _xyz_atoms(xyzname)
        elements = {sym: gto.charge(sym) for sym in {s for s, _ in atoms}}
        basis_arg = bse_json_basis(TEST / "basis_sets" / "cc-pvdz.json", elements)
        mol = gto.M(atom=atoms, basis=basis_arg, unit='Angstrom', cart=False,
                    charge=charge, spin=mult - 1, verbose=0)
        nao = mol.nao_nr()
        perm = terco_to_pyscf_perm_sph(atoms, basis_arg)
        inv_perm = np.argsort(perm)

        nspin = 1 if mult == 1 else 2
        if nspin == 1:
            def make_mf():
                return scf.RHF(mol)
            mf, grad_norm = converge_tight(make_mf, None, f"full_grad {name}")
            dm0_terco = None
        else:
            # Same symmetry-breaking recipe as gen_rks's own OH cases.
            ovlp = mol.intor('int1e_ovlp')
            dip = mol.intor('int1e_r')
            hcore_pert = scf.hf.get_hcore(mol) + 0.05 * dip[1]
            ea, ca = scipy.linalg.eigh(hcore_pert, ovlp)
            nelec = sum(gto.charge(sym) for sym, _ in atoms) - charge
            na = (nelec + (mult - 1)) // 2
            nb = nelec - na
            dm0_a = ca[:, :na] @ ca[:, :na].T
            dm0_b = ca[:, :nb] @ ca[:, :nb].T

            def make_mf():
                return scf.UHF(mol)
            mf, grad_norm = converge_tight(make_mf, (dm0_a, dm0_b), f"full_grad {name}")
            dm0_terco = np.stack([dm0_a[np.ix_(inv_perm, inv_perm)],
                                  dm0_b[np.ix_(inv_perm, inv_perm)]], axis=2)

        g = mf.nuc_grad_method()
        g.verbose = 0
        grad_py = g.kernel()   # (nat, 3), Eh/Bohr, terco's own atom order

        w.str_(f"NAME_{name}", name)
        w.int_(f"NSPIN_{name}", nspin)
        w.int_(f"NAO_{name}", nao)
        w.int_(f"NAT_{name}", len(atoms))
        if dm0_terco is not None:
            w.matrix3(f"DM0_{name}", dm0_terco)
        w.real(f"E_{name}", mf.e_tot)
        w.real(f"GRAD_NORM_{name}", grad_norm)
        w.matrix(f"GRAD_{name}", grad_py.T)   # (3, nat), matching terco's layout
        print(f"  full_grad {name}: E(pyscf,converged) {mf.e_tot:18.12f}  |g| {grad_norm:.2e}")
    w.end()
    w.write()


def gen_xc_grad(case, xyz, json_rel, functional, grid_level=3, grid_response=False):
    """The FULL KS-DFT nuclear gradient against pyscf's own grad.RKS/
    grad.UKS, at TERCO'S OWN converged orbitals and TERCO'S OWN grid
    (check_xc_grad_vs_pyscf.F90's own header has the full reasoning: no
    SCF difference, no quadrature difference, so what differs is this
    build's formula -- and the grid-response term it omits by default,
    which pyscf's grad_response=False (the default, left unset here)
    omits too).

    mo_coeff is TERCO'S cart=True AO order and normalisation; converted to
    pyscf's order (terco_to_pyscf_perm) and rescaled (ao_norm_from_overlap)
    before mf.mo_coeff/mo_occ/mo_energy are set directly -- mf.kernel() is
    NEVER called, so pyscf's own SCF never runs for this case.

    GRID_RESPONSE=TRUE is a DIFFERENT comparison, not this function's
    default path with one flag flipped: `g.grid_response=True` needs
    pyscf's OWN `grids_response_cc`, which walks `grids.atom_grids` --
    populated only by `Grids.build()`, not by overriding `coords`/
    `weights` after the fact the way the grid_response=False path below
    does with terco's dumped grid. So this path builds a plain pyscf
    `Grids(mol)` at the SAME level and lets it run its own defaults
    (Treutler-Ahlrichs radial, Bragg radii, Treutler size adjustment,
    Becke's three-iteration cutoff, NWChem pruning -- trc_dft_grid.F90's
    own defaults, confirmed to match pyscf's one by one against
    pyscf/dft/gen_grid.py and pyscf/dft/radi.py, not assumed) rather than
    terco's literal points: the comparison is then at terco's own
    orbitals, on a grid that SHOULD be the partition-for-partition
    equivalent of terco's (same scheme, same radii, same level), but not
    provably the identical list of floating-point points in the identical
    order -- see check_xc_grad_vs_pyscf.F90's own header on why that is
    the honest description of what this checks."""
    PYSCF = {'svwn': 'LDA_X,LDA_C_VWN', 'pbe': 'PBE', 'b3lyp': 'HYB_GGA_XC_B3LYP'}

    exe = f"check_xc_grad_vs_pyscf_gridresp_{case}" if grid_response else f"check_xc_grad_vs_pyscf_{case}"
    dump = TEST / f"xc_grad_probe_gridresp_{case}.bin" if grid_response else TEST / f"xc_grad_probe_{case}.bin"
    gresp_arg = [1] if grid_response else []
    run_bin(exe, [xyz, json_rel, functional, grid_level, dump.name, case] + gresp_arg)
    r = BinReader(dump)
    nat, nao, nspin, nalpha, nbeta = r.i8(), r.i8(), r.i8(), r.i8(), r.i8()
    at_z = r.i8(nat)
    rcoord = r.f8(3 * nat).reshape(3, nat, order='F').T
    npts = r.i8()
    pts = r.f8(3 * npts).reshape(3, npts, order='F').T
    wts = r.f8(npts)
    cmo = r.f8(nao * nao * nspin).reshape(nao, nao, nspin, order='F')
    eps = r.f8(nao * nspin).reshape(nao, nspin, order='F')
    grad_terco = r.f8(3 * nat).reshape(3, nat, order='F').T
    s_terco = r.mat(nao)
    r.close()

    elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
    basis_arg = bse_json_basis(TEST / json_rel, elements)
    atoms = [(_SYMBOL[int(z)], tuple(rcoord[i])) for i, z in enumerate(at_z)]
    charge = int(round(at_z.sum())) - (nalpha + nbeta)
    spin = nalpha - nbeta
    mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True,
                charge=charge, spin=spin, verbose=0)
    if mol.nao_nr() != nao:
        raise RuntimeError(f"gen_xc_grad {case}: nao mismatch terco {nao} pyscf {mol.nao_nr()}")

    perm = terco_to_pyscf_perm(atoms, basis_arg)
    inv_perm = np.argsort(perm)
    s_terco_p = s_terco[np.ix_(inv_perm, inv_perm)]
    s_pyscf = mol.intor('int1e_ovlp')
    n, s_err = ao_norm_from_overlap(s_terco_p, s_pyscf)
    if s_err >= 1e-9:
        raise RuntimeError(f"gen_xc_grad {case}: overlap check failed, {s_err:.3e}")

    if nspin == 1:
        c_p = n[:, None] * cmo[inv_perm, :, 0]
        occ = np.zeros(nao)
        occ[:nalpha] = 2.0
        mf = dft.RKS(mol)
        mf.xc = PYSCF[functional]
        if not grid_response:
            mf.grids.coords = np.ascontiguousarray(pts)
            mf.grids.weights = wts
            mf.grids.non0tab = None
        else:
            mf.grids.level = grid_level
        mf.mo_coeff = c_p
        mf.mo_occ = occ
        mf.mo_energy = eps[:, 0]
        mf.converged = True
    else:
        c_a = n[:, None] * cmo[inv_perm, :, 0]
        c_b = n[:, None] * cmo[inv_perm, :, 1]
        occ_a = np.zeros(nao)
        occ_a[:nalpha] = 1.0
        occ_b = np.zeros(nao)
        occ_b[:nbeta] = 1.0
        mf = dft.UKS(mol)
        mf.xc = PYSCF[functional]
        if not grid_response:
            mf.grids.coords = np.ascontiguousarray(pts)
            mf.grids.weights = wts
            mf.grids.non0tab = None
        else:
            mf.grids.level = grid_level
        mf.mo_coeff = (c_a, c_b)
        mf.mo_occ = (occ_a, occ_b)
        mf.mo_energy = (eps[:, 0], eps[:, 1])
        mf.converged = True

    dm = mf.make_rdm1()
    mf.e_tot = mf.energy_tot(dm)
    g = mf.nuc_grad_method()
    g.grid_response = grid_response
    #: False is the default, stated explicitly -- see this function's own
    #: docstring and trc_xc_grad.F90's header on the grid-response term
    #: both this build and pyscf's default omit. True is this function's
    #: OTHER path (see the docstring on why `mf.grids` is built fresh
    #: rather than overridden when this is set).
    grad_py = g.kernel()

    ref_tag = f"xc_grad_gridresp_{case}" if grid_response else f"xc_grad_{case}"
    w = RefWriter(f"{ref_tag}.ref")
    w.header("check_xc_grad_vs_pyscf_" + case,
              [ROOT / "test" / "check_xc_grad_vs_pyscf.F90", TEST / json_rel, TEST / xyz])
    w.int_("NAT", nat)
    w.int_("NAO", nao)
    w.matrix("GRAD", grad_py.T)
    w.end()
    w.write()
    worst = np.abs(grad_terco - grad_py).max()
    print(f"  {ref_tag}: worst|terco-pyscf| {worst:.3e} (at generation time)")
    dump.unlink()


def gen_rks():
    """scf_rks.F90's 12 cases (water RHF/SVWN/PBE/BLYP/B3LYP/wB97X-V/
    wB97X/wB97 restricted; OH radical UHF/PBE/wB97X-V/wB97X unrestricted),
    6-31G uncontracted to primitives. Commits pyscf's own converged energy
    AND its own converged density (in terco's AO order) -- the at-density
    control is then: terco evaluates ITS energy at that committed density
    (trc_scf_run with max_iter=1 and dguess=this density; see
    src/scf/trc_scf_driver.F90's own comment '# --- energy at this density ---',
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
    oh_cascade_dm = None
    #: OH's converged (dm0_a, dm0_b) in PYSCF's own AO order, carried from
    #: one open-shell case to the next (hf -> pbe -> wb97x-v -> wb97x):
    #: the symmetry-broken core guess below is strong enough to pick a
    #: DEFINITE orbital occupation for the first (hf) case, but handing
    #: every functional that SAME crude guess independently is not --
    #: wb97x from it was measured converging to a qualitatively different
    #: state (|dE| 0.149 Hartree against its own live run, not a rotation
    #: of the same one). Chaining from hf's own converged, already-
    #: resolved density keeps every functional in the same orbital
    #: character hf settled on, which a real molecule's functional-to-
    #: functional orbital ordering does not change.
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
        if nspin == 1:
            dm0 = dm_pyscf_order[:, :, 0]
            dm0_terco = None
        else:
            # OH's unpaired electron sits in one of a rotationally
            # degenerate pair of pi orbitals (the O-H axis is the only
            # distinguished direction) -- a genuine zero-curvature
            # direction of the energy, not noise in it, so EITHER
            # orientation is an equally valid converged answer and the
            # energy (and the at-density control, which only evaluates
            # the energy functional) do not care which one a run lands
            # on. The density DOES, and starting from pyscf's own default
            # guess, or even from terco's own (otherwise deterministic)
            # density, was measured landing in different rotations on
            # different invocations of this script -- not run-to-run
            # across machines, but regenerate.py run twice in the SAME
            # environment, which can only mean numpy/BLAS's own summation
            # order (thread scheduling) nudging a direction the energy
            # surface has no restoring force along. Fixed by building the
            # initial guess from the core Hamiltonian plus a SUBSTANTIAL
            # (not infinitesimal) symmetry-breaking field along one fixed
            # axis (y), diagonalized once with no SCF history behind it:
            # that guess alone picks a definite orientation strongly
            # enough that the unperturbed SCF which follows converges
            # near it rather than drifting, which an exactly-degenerate
            # starting point does not. Committed as DM0_i below so
            # test/scf_rks.F90 starts terco's OWN SCF from the identical
            # density -- both sides then converge toward the SAME
            # orientation, not merely two that happen to share an energy.
            if oh_cascade_dm is None:
                ovlp = mol.intor('int1e_ovlp')
                dip = mol.intor('int1e_r')
                hcore_pert = scf.hf.get_hcore(mol) + 0.05*dip[1]
                ea, ca = scipy.linalg.eigh(hcore_pert, ovlp)
                dm0_a = ca[:, :na] @ ca[:, :na].T
                dm0_b = ca[:, :nb] @ ca[:, :nb].T
            else:
                dm0_a, dm0_b = oh_cascade_dm
            dm0 = (dm0_a, dm0_b)
            dm0_terco = np.stack([dm0_a[np.ix_(inv_perm, inv_perm)],
                                  dm0_b[np.ix_(inv_perm, inv_perm)]], axis=2)
        mf, grad_norm = converge_tight(make_mf, dm0, f"rks case {ic} ({label})")
        e_p = mf.e_tot

        dm_conv = mf.make_rdm1()
        if nspin == 1:
            dm_conv_terco = dm_conv[np.ix_(inv_perm, inv_perm)][:, :, None]
        else:
            dm_conv_terco = np.stack([dm_conv[0][np.ix_(inv_perm, inv_perm)],
                                      dm_conv[1][np.ix_(inv_perm, inv_perm)]], axis=2)
            oh_cascade_dm = (dm_conv[0], dm_conv[1])

        w.str_(f"NAME_{ic}", name or "hf")
        w.int_(f"NSPIN_{ic}", nspin)
        w.int_(f"NAO_{ic}", n)
        # DM0_i BEFORE E_i/DM_i: test/scf_rks.F90 reads this record in two
        # passes (DM0, if present, to seed its OWN live SCF before running
        # it; the rest afterwards for the at-density control), and the
        # plain-text format's reader is strictly sequential by file
        # position, so the write order here is the read order there.
        if dm0_terco is not None:
            w.matrix3(f"DM0_{ic}", dm0_terco)
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

        # These 4 are deliberately pathological (biradicaloid N2, O3,
        # doubly-diffuse w5, an anion), and pyscf's own default guess
        # ('minao') measured right at (or just past) one gate or another
        # on more than one of them, run to run: n2_stretched's |g| sat on
        # both sides of converge_tight's 1e-8 gate across repeated runs
        # (1.539e-08, 3.477e-09, 4.397e-08, ...) from the SAME guess
        # algorithm -- not basis-conditioning noise in the gate itself,
        # but a flat enough gradient near this biradicaloid minimum that
        # minao's exact starting point decides which side of 1e-8 the
        # run lands on. 'atom' (pyscf's free-atom-density guess, not
        # minao) reaches the SAME energy on every case checked and
        # clears the gradient gate with comfortable margin every time
        # measured (n2_stretched 3.5-8.0e-09, o3 4.2e-09) -- a
        # deterministic, well-motivated guess, per Jorge's rule, not a
        # loosened gate.
        def make_mf():
            m = scf.UHF(mol) if nspin == 2 else scf.RHF(mol)
            m.init_guess = "atom"
            return m

        label = f"scf_robust {xyzname}/{basname}"
        if xyzname == "w5":
            # w5/d-aug-cc-pvdz's overlap matrix has condition number 5.3e7
            # (Cartesian, 305 AOs; smallest eigenvalue 6.5e-7 -- above
            # terco's own 1e-8 linear-dependence threshold, so terco and
            # pyscf diagonalise the identical, non-dropped basis; this is
            # not a linear-dependence bug on either side). That
            # conditioning puts a rounding floor under BOTH of
            # converge_tight's own checks on this one case: |g| measured
            # 6.675e-09 to 1.171e-08 and the self-check measured 6.9e-12
            # to 1.2e-11 across separate runs and conv_tol choices (1e-12,
            # 5e-13) with the SAME 'atom' guess -- i.e. every regeneration
            # of this one case is close enough to both gates (1e-8, 1e-11)
            # that which side it lands on is a coin toss, not a sign that
            # the SCF itself is under-converged (get_hcore is bit-
            # identical between a fresh mf and the converged one; get_veff
            # 's own trace against the converged density differs by only
            # ~1.7e-13). Tightening conv_tol does not help (checked: 5e-13
            # measured the SAME ~7e-12 self-check, not a smaller one), so
            # this is a documented exception to BOTH gates for this one
            # case only -- the other 3 keep 1e-8 / 1e-11 -- same spirit as
            # gen_scf_df's formamide/B3LYP exception a few hundred lines
            # up. This gates only PYSCF'S OWN reference quality for this
            # one hard case; it does not loosen terco's own scf_robust.F90
            # tolerance, which this branch leaves unchanged.
            mf, grad_norm = converge_tight(make_mf, None, label, conv_tol_grad=5.0e-8, self_check_tol=2.0e-11)
        else:
            mf, grad_norm = converge_tight(make_mf, None, label)
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


def _xyz_atoms(xyzname):
    """Plain (symbol, (x, y, z)) list off a raw .xyz file, Angstrom, no
    Bohr conversion -- the same read scf_robust's own gen_scf_robust does,
    for `gto.M(..., unit='Angstrom')`."""
    atoms = []
    with open(TEST / (xyzname + ".xyz")) as g:
        nat = int(g.readline())
        g.readline()
        for _ in range(nat):
            parts = g.readline().split()
            atoms.append((parts[0], (float(parts[1]), float(parts[2]), float(parts[3]))))
    return atoms


def _gen_scf_df_one(r, name, xyzname, basname):
    """One of test/scf_df.F90's 12 cases, read off the shared BinReader
    `r`: nao/natm/nalpha/nbeta (a sanity check only), the molecule and
    basis names it used (also passed in, read from the probe for the
    SAME check), the method name and terco's own converged energy.

    Rebuilds the orbital and jkfit (def2-universal-jkfit) bases itself,
    via `bse_json_basis` straight off the BSE JSON -- NOT from terco's own
    dumped sh_c/sh_e, which fold in a primitive normalisation pyscf's
    custom-basis-dict path would apply a second time (see scf_df.F90's
    own header for how this was found: every multi-primitive contraction
    came out wrong fed that way, matching bse_json_basis exactly once
    reconstructed from the JSON instead)."""
    PYSCF = {'pbe': 'PBE', 'b3lyp': 'HYB_GGA_XC_B3LYP'}
    n, natm, na, nb = (r.i8() for _ in range(4))
    mol_ref = r.str_(16)
    bas_ref = r.str_(16)
    method = r.str_(8)
    e_terco = r.f8(1)[0]
    if mol_ref != xyzname or bas_ref != basname:
        raise RuntimeError(f"gen_scf_df {name}: probe order mismatch "
                            f"({mol_ref}/{bas_ref} vs {xyzname}/{basname})")

    atoms = _xyz_atoms(xyzname)
    if len(atoms) != natm:
        raise RuntimeError(f"gen_scf_df {name}: natm mismatch terco {natm} xyz {len(atoms)}")
    elements = {sym: {'H': 1, 'C': 6, 'N': 7, 'O': 8}[sym] for sym, _ in atoms}
    bas = bse_json_basis(TEST / "basis_sets" / (basname + ".json"), elements)
    auxbas = bse_json_basis(TEST / "basis_sets" / "def2-universal-jkfit.json", elements)
    mol = gto.M(atom=atoms, basis=bas, cart=True, unit='Angstrom', verbose=0)
    if mol.nao_nr() != n:
        raise RuntimeError(f"gen_scf_df {name}: nao mismatch terco {n} pyscf {mol.nao_nr()}")

    def make_mf(grid_level=3):
        if method:
            m = dft.RKS(mol)
            m.xc = PYSCF[method]
            m.grids.level = grid_level
        else:
            m = scf.RHF(mol)
        return m.density_fit(auxbasis=auxbas)

    label = f"scf_df {name} ({method or 'hf'})"
    try:
        mf, grad_norm = converge_tight(make_mf, None, label)
    except ConvergenceError:
        # formamide/cc-pVDZ/B3LYP measured: converged=True but |g| sits at
        # ~1.1-1.3e-08 regardless of the XC grid level (3, 5, 6, 9 all
        # measured) or a .newton() retry -- a quadrature/DF floor in the
        # gradient on this one case, not a missed stationary point. One
        # documented exception, same spirit as rks.ref's own looser bound
        # for a hard open-shell DFT case (test/pyscf_refs/README.md): this
        # case alone, at 5e-8 (5x tighter than pyscf's own sqrt(conv_tol)
        # default of ~1e-6).
        mf, grad_norm = converge_tight(make_mf, None, label + " (loosened gate)", conv_tol_grad=5.0e-8)
    e_p = mf.e_tot
    print(f"  {name:28s} E(terco) {e_terco:20.12f}  E(pyscf) {e_p:20.12f}  "
          f"|dE| {abs(e_terco - e_p):.2e}  |g| {grad_norm:.2e}")
    return e_p, mf.conv_tol, mf.conv_tol_grad, grad_norm


def gen_scf_df():
    """test/scf_df.F90's 12 cases (water/formamide x cc-pVDZ/6-31G* x
    HF/PBE/B3LYP, all --df) against pyscf's own RHF/RKS.density_fit on the
    identical auxiliary basis (def2-universal-jkfit, mol.cart=True)."""
    run_bin("scf_df", [])
    dump = TEST / "scf_df_probe.bin"
    r = BinReader(dump)
    ncase = r.i8()
    w = RefWriter("scf_df.ref")
    w.header("scf_df", [ROOT / "test" / "scf_df.F90",
                        TEST / "basis_sets" / "def2-universal-jkfit.json",
                        TEST / "basis_sets" / "cc-pvdz.json", TEST / "basis_sets" / "6-31gs.json",
                        TEST / "water.xyz", TEST / "formamide.xyz"])
    w.int_("NCASE", ncase)
    gmol = ["water", "water", "formamide", "formamide"]
    gbasis = ["cc-pvdz", "6-31gs", "cc-pvdz", "6-31gs"]
    methods = ["hf", "pbe", "b3lyp"]
    ic = 0
    for ig in range(4):
        for im in range(3):
            name = f"{gmol[ig]}_{gbasis[ig]}_{methods[im]}"
            e_p, ctol, ctolg, gnorm = _gen_scf_df_one(r, name, gmol[ig], gbasis[ig])
            w.str_(f"NAME_{ic}", name)
            w.real(f"E_{ic}", e_p)
            w.real(f"CONV_TOL_{ic}", ctol)
            w.real(f"CONV_TOL_GRAD_{ic}", ctolg)
            w.real(f"GRAD_NORM_{ic}", gnorm)
            ic += 1
    r.close()
    w.end()
    w.write()
    dump.unlink()


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
    lines = ["# kernel family polar npts nin nout", "# libxc %s via pyscf" % libxc.libxc_version()]
    lines += ["# " + h for h in ref_header_lines([ROOT / "tools" / "xc_kernel_discovery.py"])]
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


# =========================================================================
# The polarisable continuum (src/solvation/trc_pcm.F90 and the --pcm kernels).
#
# Every case uses pyscf.solvent.pcm's own machinery -- gen_surface, get_D_S,
# the PCM object's K/R, _get_v/_get_vmat (int3c2e against
# fakemol_for_charges, pyscf's own surface-integral path) -- so what is
# committed is pyscf's continuum, not a re-derivation of it.
# =========================================================================

PCM_EPS = 78.3553          #: pyscf's PCM default (water); also the tests'
PCM_ANG = 1.8897261254578281
#: trc_test_basis.F90's ANG: the cavity cases read water.xyz themselves
#: (there is no binary to run for them), so they convert Angstrom to bohr
#: with the SAME constant read_xyz does, and hand pyscf bohr.


def _xyz_bohr(xyzname):
    atoms = []
    with open(TEST / xyzname) as g:
        nat = int(g.readline())
        g.readline()
        for _ in range(nat):
            parts = g.readline().split()
            atoms.append((parts[0], tuple(float(c) * PCM_ANG for c in parts[1:4])))
    return atoms


def gen_pcm_cavity():
    """check_pcm_cavity.F90: water's cavity at 302 points per sphere (points,
    weights, exponents, switching values, areas, normals; S x and D x for a
    fixed x; both models' symmetrised charges and q.v/2 for a fixed v), and
    S and D element by element at 110 points per sphere, small enough to
    commit whole. The fixed x and v are the closed forms check_pcm_cavity.F90
    evaluates on its own points (xfix/vfix there) -- functions of position,
    so pyscf's point order needs no matching here; the Fortran side matches
    points by position before comparing."""
    from pyscf.solvent import pcm
    atoms = _xyz_bohr("water.xyz")
    mol = gto.M(atom=atoms, basis='sto-3g', unit='Bohr', verbose=0)
    rad = pcm.modified_Bondi * 1.2

    def vfix(c):
        return (1.0 / np.linalg.norm(c - [0.1, 0.2, 0.3], axis=1)
                + 0.5 / np.linalg.norm(c - [-0.5, 0.4, -0.2], axis=1))

    def xfix(c):
        return np.cos(1.3 * c[:, 0]) + np.sin(0.7 * c[:, 1]) + 0.2 * c[:, 2]

    surf = pcm.gen_surface(mol, ng=302, rad=rad)
    c = surf['grid_coords']
    D, S = pcm.get_D_S(surf, with_S=True, with_D=True)
    x = xfix(c)
    v = vfix(c)

    w = RefWriter("pcm_cavity_water.ref")
    w.header("check_pcm_cavity", [ROOT / "test" / "check_pcm_cavity.F90", TEST / "water.xyz"])
    w.int_("NPTS", c.shape[0])
    w.matrix("COORDS", c.T)
    w.vector("WEIGHTS", surf['weights'])
    w.vector("ZETA", surf['charge_exp'])
    w.vector("SWITCH", surf['switch_fun'])
    w.vector("AREA", surf['area'])
    w.matrix("NORMALS", surf['norm_vec'].T)
    w.vector("SX", S @ x)
    w.vector("DX", D @ x)
    for method, tag in (("C-PCM", "CPCM"), ("IEF-PCM", "IEFPCM")):
        sol = pcm.PCM(mol)
        sol.method = method
        sol.eps = PCM_EPS
        sol.lebedev_order = 29
        sol.vdw_scale = 1.2
        sol.build()
        assert np.allclose(sol.surface['grid_coords'], c, atol=0, rtol=0)
        K, R = sol._intermediates['K'], sol._intermediates['R']
        q = np.linalg.solve(K, R @ v)
        qt = R.T @ np.linalg.solve(K.T, v)
        qs = 0.5 * (q + qt)
        w.vector(f"Q_{tag}", qs)
        w.real(f"E_{tag}", 0.5 * qs @ v)
        print(f"  pcm_cavity: {method:8s} {c.shape[0]} points, q.v/2 = {0.5 * qs @ v:.12f}")

    surf = pcm.gen_surface(mol, ng=110, rad=rad)
    D, S = pcm.get_D_S(surf, with_S=True, with_D=True)
    w.int_("NPTS110", surf['grid_coords'].shape[0])
    w.matrix("COORDS110", surf['grid_coords'].T)
    w.matrix("S110", S)
    w.matrix("D110", D)
    w.end()
    w.write()


def _pcm_v_nj(mol, coords, zeta):
    """(nao, nao, npts) <i|erf(zeta_k|r - C_k|)/|r - C_k||j>, pyscf.solvent.pcm's
    own way (PCM._get_v/_get_vmat: int3c2e against fakemol_for_charges)."""
    from pyscf.gto import mole
    fakemol = mole.fakemol_for_charges(coords, expnt=zeta**2)
    fakemol.cart = mol.cart
    return df.incore.aux_e2(mol, fakemol, intor='int3c2e', aosym='s1')


def gen_pcm_ints(case, xyz, json_rel):
    """check_pcm_ints.F90: the PCM surface-integral kernels on fixed inputs
    (its own cavity points plus eight off-cavity ones, a fixed density,
    fixed charges). Committed: v_k = D . I_k and F = sum_k q_k I_k from
    pyscf's int3c2e/fakemol path, after checking that path against
    int1e_grids under mol.with_range_coulomb(zeta), one zeta group at a
    time -- the two pyscf routes to the same erf(zeta r)/r integral."""
    dump = TEST / f"pcm_ints_probe_{case}.bin"
    run_bin("check_pcm_ints_" + case, [xyz, json_rel, dump.name, case])
    if not dump.exists():
        print(f"  pcm_ints {case}: binary self-skipped (LMAX too low for this build); not regenerating")
        return
    r = BinReader(dump)
    nat, nao, npts = r.i8(), r.i8(), r.i8()
    at_z = r.i8(nat)
    rcoord = r.f8(3 * nat).reshape(3, nat, order='F').T
    pts = r.f8(3 * npts).reshape(3, npts, order='F').T
    zeta = r.f8(npts)
    q = r.f8(npts)
    d_terco = r.mat(nao)
    s_terco = r.mat(nao)
    r.close()

    elements = {_SYMBOL[int(z)]: int(z) for z in at_z}
    basis_arg = bse_json_basis(TEST / json_rel, elements)
    atoms = [(_SYMBOL[int(z)], tuple(rcoord[i])) for i, z in enumerate(at_z)]
    mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=True, verbose=0)
    if mol.nao_nr() != nao:
        raise RuntimeError(f"gen_pcm_ints {case}: nao mismatch terco {nao} pyscf {mol.nao_nr()}")
    perm = terco_to_pyscf_perm(atoms, basis_arg)
    inv_perm = np.argsort(perm)
    s_terco_p = s_terco[np.ix_(inv_perm, inv_perm)]
    d_terco_p = d_terco[np.ix_(inv_perm, inv_perm)]
    n, s_err = ao_norm_from_overlap(s_terco_p, mol.intor('int1e_ovlp'))
    if s_err >= 1e-10:
        raise RuntimeError(f"gen_pcm_ints {case}: overlap check failed, {s_err:.3e}")
    d_pyscf = n[:, None] * d_terco_p * n[None, :]

    v_nj = _pcm_v_nj(mol, pts, zeta)
    worst = 0.0
    for z in np.unique(zeta):
        idx = np.where(zeta == z)[0]
        with mol.with_range_coulomb(z):
            g = mol.intor('int1e_grids', grids=pts[idx])        # (ng, nao, nao)
        worst = max(worst, np.abs(np.transpose(g, (1, 2, 0)) - v_nj[:, :, idx]).max())
    if worst > 1e-12:
        raise RuntimeError(f"gen_pcm_ints {case}: int3c2e/fakemol and int1e_grids "
                           f"disagree by {worst:.3e}")
    print(f"  pcm_ints {case}: int3c2e/fakemol vs int1e_grids(range_coulomb): {worst:.2e}")
    v = np.einsum('ijL,ij->L', v_nj, d_pyscf)
    f_pyscf = np.einsum('ijL,L->ij', v_nj, q)
    f_terco = (f_pyscf * n[:, None] * n[None, :])[np.ix_(perm, perm)]

    w = RefWriter(f"pcm_ints_{case}.ref")
    w.header("check_pcm_ints_" + case,
             [ROOT / "test" / "check_pcm_ints.F90", TEST / json_rel, TEST / xyz])
    w.int_("NAO", nao)
    w.int_("NPTS", npts)
    w.matrix("S", s_terco)
    w.vector("V", v)
    w.matrix("F", f_terco)
    w.end()
    w.write()
    dump.unlink()


PCM_CONV_TOL = 1e-12
PCM_CONV_TOL_GRAD = 1e-8


def _pcm_converge(mf, label, dm0=None):
    """Converge one PCM SCF reference tightly or refuse it: conv_tol 1e-12
    AND conv_tol_grad 1e-8, a second-order (.newton()) retry when the
    first-order solver does not get there, then two self-checks -- the
    orbital-gradient norm actually under the gate, and energy_tot at the
    converged density reproducing e_tot (the at-density control compares
    against the former and the converged energy against the latter, so
    they must be the same number). Returns (mf, grad_norm)."""
    mf.conv_tol = PCM_CONV_TOL
    mf.conv_tol_grad = PCM_CONV_TOL_GRAD
    mf.max_cycle = 300
    mf.verbose = 0
    mf.kernel(dm0=dm0)
    if not mf.converged:
        mo, occ = mf.mo_coeff, mf.mo_occ
        mf = mf.newton()
        mf.conv_tol = PCM_CONV_TOL
        mf.conv_tol_grad = PCM_CONV_TOL_GRAD
        mf.max_cycle = 300
        mf.kernel(mo, occ)
    require_converged(mf, label)

    def grad_norm():
        # The Fock matrix handed in explicitly: get_grad's own fallback
        # builds it from get_hcore + get_veff, which leaves out the
        # solvent's potential (SCFWithSolvent adds it in get_fock).
        fock = mf.get_fock(dm=mf.make_rdm1())
        return float(np.linalg.norm(mf.get_grad(mf.mo_coeff, mf.mo_occ, fock)))

    # pyscf's own test reads the gradient of the step BEFORE the last one,
    # so the final density can sit just above conv_tol_grad with the flag
    # set. Then keep going, from there, at a tenfold tighter gate -- never
    # accept the larger number.
    gnorm = grad_norm()
    for _ in range(3):
        if gnorm <= PCM_CONV_TOL_GRAD:
            break
        mf.conv_tol_grad = 0.1 * PCM_CONV_TOL_GRAD
        mf.kernel(dm0=mf.make_rdm1())
        require_converged(mf, label)
        gnorm = grad_norm()
    if gnorm > PCM_CONV_TOL_GRAD:
        raise ConvergenceError(f"{label}: |g| = {gnorm:.2e} > {PCM_CONV_TOL_GRAD} after converging")
    de = abs(mf.energy_tot(mf.make_rdm1()) - mf.e_tot)
    if de > 1e-10:
        raise ConvergenceError(f"{label}: energy_tot(dm) differs from e_tot by {de:.2e}")
    return mf, gnorm


def gen_pcm_scf(set_name):
    """scf_pcm.F90's cases: pyscf's own RHF/UHF/RKS/UKS with .PCM() (same
    model, dielectric, radii scale and Lebedev order), on terco's own DFT
    grid (and VV10 grid), converged by _pcm_converge (conv_tol 1e-12 AND
    conv_tol_grad 1e-8, refused if not). Commits the converged total energy and E_diel
    (scf_summary['e_solvent'], the one that went into e_tot), with each
    case's orbital-gradient norm in a comment, and, for every
    Cartesian case, pyscf's converged density in terco's AO order and
    normalisation with energy_tot(dm) and E_diel at that density -- the
    at-density control (see scf_rks.F90/gen_rks for the pattern)."""
    from pyscf.dft import gen_grid
    PYSCF_XC = {'b3lyp': 'HYB_GGA_XC_B3LYP', 'wb97x-v': 'wb97x-v', 'pbe': 'PBE'}
    PYSCF_PCM = {'cpcm': 'C-PCM', 'iefpcm': 'IEF-PCM'}
    LEB_ORDER = {npts: order for order, npts in gen_grid.LEBEDEV_ORDER.items()}

    run_bin("scf_pcm", [set_name])
    dump = TEST / f"pcm_scf_probe_{set_name}.bin"
    if not dump.exists():
        print(f"  pcm_scf {set_name}: binary self-skipped (LMAX too low for this build); not regenerating")
        return
    r = BinReader(dump)
    ncase = r.i8()
    w = RefWriter(f"pcm_scf_{set_name}.ref")
    w.header("scf_pcm " + set_name,
             [ROOT / "test" / "scf_pcm.F90", TEST / "water.xyz", TEST / "formamide.xyz", TEST / "oh.xyz",
              TEST / "nh2.xyz",
              TEST / "basis_sets" / "6-31gs.json", TEST / "basis_sets" / "cc-pvdz.json",
              TEST / "basis_sets" / "cc-pvtz.json"])
    w.comment(f"CONV_TOL {PCM_CONV_TOL:g}")
    w.comment(f"CONV_TOL_GRAD {PCM_CONV_TOL_GRAD:g}")
    w.int_("NCASE", ncase)
    for ic in range(ncase):
        nao_c, nao_scf, natm, nspin, na, nb, sph = (r.i8() for _ in range(7))
        xyzname, basname, fn, model = r.str_(16), r.str_(32), r.str_(8), r.str_(8)
        eps = r.f8(1)[0]
        leb = r.i8()
        rscale = r.f8(1)[0]
        zint = r.i8(natm)
        rr = r.f8(3 * natm).reshape(3, natm, order='F').T
        npts = r.i8()
        if npts:
            pts = r.f8(3 * npts).reshape(3, npts, order='F').T
            wts = r.f8(npts)
        npts_nlc = r.i8()
        if npts_nlc:
            pts_nlc = r.f8(3 * npts_nlc).reshape(3, npts_nlc, order='F').T
            wts_nlc = r.f8(npts_nlc)
        s_terco = r.mat(nao_c)
        e_t, ep_t = r.f8(2)
        dm_t_conv = r.f8(nao_scf * nao_scf * nspin).reshape(nao_scf, nao_scf, nspin, order='F')

        elements = {_SYMBOL[int(z)]: int(z) for z in zint}
        basis_arg = bse_json_basis(TEST / "basis_sets" / (basname + ".json"), elements)
        atoms = [(_SYMBOL[int(z)], tuple(rr[i])) for i, z in enumerate(zint)]
        charge = int(sum(zint)) - (na + nb)
        mol = gto.M(atom=atoms, basis=basis_arg, unit='Bohr', cart=not sph, charge=charge,
                    spin=na - nb, verbose=0)
        if mol.nao_nr() != nao_scf:
            raise RuntimeError(f"gen_pcm_scf case {ic}: nao mismatch terco {nao_scf} pyscf {mol.nao_nr()}")

        if fn:
            mf = dft.RKS(mol) if nspin == 1 else dft.UKS(mol)
            mf.xc = PYSCF_XC[fn]
        else:
            mf = scf.RHF(mol) if nspin == 1 else scf.UHF(mol)
        mf = mf.PCM()
        mf.with_solvent.method = PYSCF_PCM[model]
        mf.with_solvent.eps = eps
        mf.with_solvent.lebedev_order = LEB_ORDER[leb]
        mf.with_solvent.vdw_scale = rscale
        if fn:
            mf.grids.coords = np.ascontiguousarray(pts)
            mf.grids.weights = wts
            mf.grids.non0tab = None
            if npts_nlc:
                mf.nlcgrids.coords = np.ascontiguousarray(pts_nlc)
                mf.nlcgrids.weights = wts_nlc
                mf.nlcgrids.non0tab = None
        label = f"{xyzname}/{basname} {fn or 'hf'} {model} {'U' if nspin == 2 else 'R'}"
        # Started from terco's converged density, as gen_rks starts its
        # cases: the OH radical's degenerate pi pair leaves UHF stationary
        # points 3e-8 Eh apart in a Lebedev cavity (4e-7 for B3LYP), and the
        # comparison is only meaningful within one basin (scf_pcm.F90's
        # header). pyscf still converges on its own to conv_tol_grad; the
        # start only picks the basin. A spherical
        # case has no AO map to pyscf's spherical convention here and keeps
        # pyscf's own guess (closed-shell water: one basin).
        dm0 = None
        if not sph:
            perm = terco_to_pyscf_perm(atoms, basis_arg)
            inv_perm = np.argsort(perm)
            n, s_err = ao_norm_from_overlap(s_terco[np.ix_(inv_perm, inv_perm)], mol.intor('int1e_ovlp'))
            if s_err >= 1e-10:
                raise RuntimeError(f"gen_pcm_scf case {ic}: overlap check failed, {s_err:.3e}")
            nn = n[:, None] * n[None, :]
            dmp = [nn * dm_t_conv[:, :, k][np.ix_(inv_perm, inv_perm)] for k in range(nspin)]
            dm0 = dmp[0] if nspin == 1 else np.stack(dmp)
        mf, gnorm = _pcm_converge(mf, f"pcm_scf {set_name} case {ic} ({label})", dm0)
        e_p = mf.e_tot
        ep_p = mf.scf_summary['e_solvent']

        w.comment(f"case {ic}: {label}, GRAD_NORM {gnorm:.3e}")
        w.int_(f"NSPIN_{ic}", nspin)
        w.int_(f"NAO_{ic}", nao_scf)
        w.real(f"E_{ic}", e_p)
        w.real(f"EPCM_{ic}", ep_p)
        extra = ""
        if not sph:
            dm = mf.make_rdm1()
            e_dm = mf.energy_tot(dm)
            ep_dm = mf.with_solvent.kernel(dm)[0]
            dms = [dm] if nspin == 1 else [dm[0], dm[1]]
            dm_t = np.stack([(d / nn)[np.ix_(perm, perm)] for d in dms], axis=2)
            w.real(f"EDM_{ic}", e_dm)
            w.real(f"EPCMDM_{ic}", ep_dm)
            w.matrix3(f"DM_{ic}", dm_t)
            extra = f"  E(dm) - E {e_dm - e_p:.1e}"
        print(f"  pcm_scf {ic:2d} {label:34s} E {e_p:18.12f} (terco {e_t - e_p:+.1e})  "
              f"E_diel {ep_p:.12f} (terco {ep_t - ep_p:+.1e}){extra}")
    w.end()
    w.write()
    dump.unlink()


CASES = {
    "mult": gen_mult,
    "ao": gen_ao,
    "xc_energy": gen_xc_energy,
    "vv10": gen_vv10,
    "rys_grad_631gs": lambda: gen_rys_grad("631gs", "basis_sets/6-31gs.json"),
    "rys_grad_ccpvdz": lambda: gen_rys_grad("ccpvdz", "basis_sets/cc-pvdz.json"),
    "onee_grad_water_631gs": lambda: gen_1e_grad("water_631gs", "water.xyz", "basis_sets/6-31gs.json"),
    "onee_grad_water_ccpvdz": lambda: gen_1e_grad("water_ccpvdz", "water.xyz", "basis_sets/cc-pvdz.json"),
    "onee_grad_formamide_ccpvdz": lambda: gen_1e_grad(
        "formamide_ccpvdz", "formamide.xyz", "basis_sets/cc-pvdz.json"),
    "xc_grad_water_svwn_631gs": lambda: gen_xc_grad(
        "water_svwn_631gs", "water.xyz", "basis_sets/6-31gs.json", "svwn"),
    "xc_grad_water_pbe_631gs": lambda: gen_xc_grad(
        "water_pbe_631gs", "water.xyz", "basis_sets/6-31gs.json", "pbe"),
    "xc_grad_water_b3lyp_631gs": lambda: gen_xc_grad(
        "water_b3lyp_631gs", "water.xyz", "basis_sets/6-31gs.json", "b3lyp"),
    "xc_grad_water_svwn_ccpvdz": lambda: gen_xc_grad(
        "water_svwn_ccpvdz", "water.xyz", "basis_sets/cc-pvdz.json", "svwn"),
    "xc_grad_water_pbe_ccpvdz": lambda: gen_xc_grad(
        "water_pbe_ccpvdz", "water.xyz", "basis_sets/cc-pvdz.json", "pbe"),
    "xc_grad_water_b3lyp_ccpvdz": lambda: gen_xc_grad(
        "water_b3lyp_ccpvdz", "water.xyz", "basis_sets/cc-pvdz.json", "b3lyp"),
    "xc_grad_formamide_pbe_ccpvdz": lambda: gen_xc_grad(
        "formamide_pbe_ccpvdz", "formamide.xyz", "basis_sets/cc-pvdz.json", "pbe"),
    "xc_grad_formamide_b3lyp_ccpvdz": lambda: gen_xc_grad(
        "formamide_b3lyp_ccpvdz", "formamide.xyz", "basis_sets/cc-pvdz.json", "b3lyp"),
    "xc_grad_oh_pbe_631gs": lambda: gen_xc_grad(
        "oh_pbe_631gs", "oh.xyz", "basis_sets/6-31gs.json", "pbe"),
    "xc_grad_gridresp_water_svwn_631gs": lambda: gen_xc_grad(
        "water_svwn_631gs", "water.xyz", "basis_sets/6-31gs.json", "svwn", grid_response=True),
    "xc_grad_gridresp_water_pbe_631gs": lambda: gen_xc_grad(
        "water_pbe_631gs", "water.xyz", "basis_sets/6-31gs.json", "pbe", grid_response=True),
    "xc_grad_gridresp_nh2_pbe_631gs": lambda: gen_xc_grad(
        "nh2_pbe_631gs", "nh2.xyz", "basis_sets/6-31gs.json", "pbe", grid_response=True),
    "rys_fock_water_ccpvdz": lambda: gen_rys_fock("water_ccpvdz", "water.xyz", "basis_sets/cc-pvdz.json", 1e-9),
    "rys_fock_water_ccpvtz": lambda: gen_rys_fock("water_ccpvtz", "water.xyz", "basis_sets/cc-pvtz.json", 1e-10),
    "rys_fock_formamide_def2tzvp": lambda: gen_rys_fock(
        "formamide_def2tzvp", "formamide.xyz", "basis_sets/def2-tzvp.json", 1e-9),
    "erf_fock_water_631gs": lambda: gen_erf_fock("water_631gs", "water.xyz", "basis_sets/6-31gs.json"),
    "erf_fock_water_ccpvdz": lambda: gen_erf_fock("water_ccpvdz", "water.xyz", "basis_sets/cc-pvdz.json"),
    "erf_fock_formamide_ccpvdz": lambda: gen_erf_fock("formamide_ccpvdz", "formamide.xyz", "basis_sets/cc-pvdz.json"),
    "df_vs_pyscf": gen_df_vs_pyscf,
    "full_grad": gen_full_grad,
    "rks": gen_rks,
    "scf_robust": gen_scf_robust,
    "rimp2": gen_rimp2,
    "scf_df": gen_scf_df,
    "xc_kernels": gen_xc_kernels,
    "pcm_cavity": gen_pcm_cavity,
    "pcm_ints_water_631gs": lambda: gen_pcm_ints("water_631gs", "water.xyz", "basis_sets/6-31gs.json"),
    "pcm_ints_water_ccpvtz": lambda: gen_pcm_ints("water_ccpvtz", "water.xyz", "basis_sets/cc-pvtz.json"),
    "pcm_scf_main": lambda: gen_pcm_scf("main"),
    "pcm_scf_f": lambda: gen_pcm_scf("f"),
}


def main():
    global _CURRENT_CASE
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
        _CURRENT_CASE = name
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
