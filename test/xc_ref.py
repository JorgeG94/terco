#!/usr/bin/env python3
"""Record libxc's values for the translated XC kernels, through PySCF.

Writes `reference_xc.dat`: for each kernel, unpolarised and polarised, a
block of points with the inputs libxc was given and every derivative it
returned through second order. `check_xc.F90` reads the file, evaluates
terco's translated kernels at the same points, and compares.

The points are chosen to cover what a real grid produces: densities from
1e-9 to 10 and reduced gradients from zero to far past the enhancement
factor's saturation, plus spin polarisations up to 0.99 with the cross
gradient at several angles. Nothing sits at a screening threshold: the
translated kernels carry ExchCXX's thresholds and libxc its own, and where
they differ the two are entitled to disagree.

The lower end is 1e-9, not 1e-10, because libxc is the one that gives out
first: below a SPIN density of about 1e-11 its polarised GGA potentials
drift from the closed form (1e-12 relative at 5e-12, 1e-10 at 5e-13), and
its polarised second derivatives carry cancellation residue of 1e-9. The
translated kernels match the closed form to 1e-14 all the way down. See
the tolerance note in check_xc.F90.

Usage:
    python test/xc_ref.py            # writes test/reference_xc.dat
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

try:
    from pyscf.dft import libxc
except ImportError:
    print("pyscf is not importable; reference_xc.dat left as it is", file=sys.stderr)
    sys.exit(0)

# terco kernel stem -> libxc name
KERNELS = {
    "slater_exchange": ("LDA_X", "lda"),
    "vwn": ("LDA_C_VWN", "lda"),
    "vwn_rpa": ("LDA_C_VWN_RPA", "lda"),
    "pbe_x": ("GGA_X_PBE", "gga"),
    "pbe_c": ("GGA_C_PBE", "gga"),
    "b88": ("GGA_X_B88", "gga"),
    "lyp": ("GGA_C_LYP", "gga"),
    # HYB_GGA_XC_WB97X_V's libxc KERNEL is purely semilocal (GGA): the
    # omega-attenuated SR exchange hole model, evaluated analytically, plus
    # VV10's semilocal correlation partner. The exact-exchange mixing
    # (exx/exx_lr in trc_xc_functional) and VV10's nonlocal term are both
    # OUTSIDE this kernel, exactly as for pbe_x inside pbe0.
    "wb97x_v_xc": ("WB97X_V", "gga"),
}

RHO = np.logspace(-9, 1, 11)
# sigma = s * rho^(8/3): s is the square of the reduced gradient, up to scale.
S_RED = np.array([0.0, 1e-3, 0.1, 1.0, 10.0, 1e3])
ZETA = np.array([0.0, 0.3, 0.7, 0.99])
COS_AB = np.array([-0.8, 0.0, 0.6])

# wb97x_v_xc's exchange part is Iikura-Tsuneda-Yanai's short-range exchange
# hole, built from erf/exp of a = omega/(2*kF); for density rho (unpolarised)
# or 2*rho_a / 2*rho_b per spin channel (polarised, by the exchange spin-
# scaling relation), the kernel MERGEs to an asymptotic series once a >=
# 1.35 to avoid cancellation in the direct erf/exp form. The switch sits at
#   rho* = (C / 1.35)**3,  C = (9/pi)**(2/3) * omega * 3**(1/3) / 18
# with omega = 0.3 -- about 4.633e-5 (unpolarised) or 2.316e-5 per spin
# channel (polarised). Points straddling it exercise both branches AND the
# boundary itself, where the two closed forms must agree to rounding.
_C_SWITCH = (9.0 ** (2.0 / 3.0)) * (0.682784063255295503581 ** 2) * 0.3 * 1.442249570307408382321638310780109588392 / 18.0
RHO_SWITCH = (_C_SWITCH / 1.35) ** 3
RHO_A_SWITCH = RHO_SWITCH / 2.0
_SWITCH_DELTA = np.array([-1e-3, -1e-6, 0.0, 1e-6, 1e-3])


def wb97x_v_extra_points(polar):
    """Points straddling the omega-attenuation branch switch, not on the
    generic RHO/S_RED/ZETA grid (which does not land near rho* exactly)."""
    if not polar:
        pts = [(rho, s * rho ** (8.0 / 3.0))
               for rho in RHO_SWITCH * (1.0 + _SWITCH_DELTA) for s in S_RED]
        return np.array(pts)
    pts = []
    other = 1.0e-3  # a density far from its own switch (2*other << rho*)
    for ra_is_switch in (True, False):
        for d in _SWITCH_DELTA:
            ra = RHO_A_SWITCH * (1.0 + d) if ra_is_switch else other
            rb = other if ra_is_switch else RHO_A_SWITCH * (1.0 + d)
            for s in S_RED[::2]:
                saa = s * ra ** (8.0 / 3.0)
                sbb = s * rb ** (8.0 / 3.0)
                for c in COS_AB:
                    pts.append((ra, rb, saa, c * np.sqrt(saa * sbb), sbb))
    return np.array(pts)


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
    """PySCF's rho layout: (rho, dx, dy, dz) per spin for a GGA."""
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
    # grad_b at an angle to grad_a chosen so that grad_a . grad_b = sab
    gb = np.sqrt(sbb)
    with np.errstate(invalid="ignore", divide="ignore"):
        cos = np.where(a[1] * gb > 0, sab / (a[1] * gb), 0.0)
    cos = np.clip(cos, -1.0, 1.0)
    b[1] = gb * cos
    b[2] = gb * np.sqrt(1.0 - cos**2)
    return np.stack([a, b])


def main():
    out = Path(__file__).with_name("reference_xc.dat")
    lines = ["# kernel family polar npts nin nout", "# libxc %s via pyscf" % libxc.libxc_version()]
    for stem, (name, family) in KERNELS.items():
        for polar in (0, 1):
            pts = polar_points(family) if polar else unpolar_points(family)
            if stem == "wb97x_v_xc":
                pts = np.concatenate([pts, wb97x_v_extra_points(polar == 1)])
            rho = libxc_rho(family, polar, pts)
            exc, vxc, fxc, _ = libxc.eval_xc(name, rho, spin=polar, deriv=2)
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
            for p, r in zip(pts, ref):
                lines.append(" ".join(f"{x:.17e}" for x in np.concatenate([p, r])))
            print(f"{stem:16s} polar={polar} {len(pts):4d} points, {ref.shape[1]:2d} outputs")
    out.write_text("\n".join(lines) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
