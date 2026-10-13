#!/usr/bin/env python3
"""Generate a Rys roots/weights table for one root count.

WHY A TABLE AND NOT A CONSTRUCTION
----------------------------------
GAMESS builds the roots per quartet: a fixed 35-point base rule reweighted by
exp(-X u), discretised Stieltjes for the Jacobi matrix, then Golub-Welsch by QL.
That is correct and needs only double precision, but it is unusable on the
device this library targets: compiled for it, that construction alone reports
REG:255, the hardware ceiling, because its 4x35 doubles of Stieltjes scratch
consume the whole register budget before a single 2D integral exists. Add 35
exponentials and a QL whose trip count is data dependent -- warp divergence --
and the case for tabulating is overwhelming.

Interpolation costs ~130 FMAs and holds only its 2n outputs.

ACCURACY AND SIZE
-----------------
Piecewise Chebyshev in X on [0, 55], degree 13, interval width 2.5. Measured
max relative error 7.0e-15 on the roots and 9.5e-15 on the weights -- the
double precision floor. 22 intervals x 2n series x 14 coefficients.

Width 2.5 rather than 5 because of the WEIGHTS, not the roots. At width 5 the
roots already reach 6.4e-15 but the weights only manage 5.3e-10, and raising
the degree barely helps (width 5 degree 17 still gives 6.5e-14) while halving
the width fixes it outright. The weight functions have a smaller analyticity
radius than the roots, so more intervals beat a higher degree.

Above X = 55 the asymptotic form u_i -> c_i/X, w_i -> d_i/sqrt(X) is exact to
double precision, so the table stops there. GAMESS switches at the same point.

The reference values come from mpmath at 40 digits -- moments, Chebyshev
algorithm for the Jacobi matrix, symmetric eigenproblem -- NOT from a double
precision construction, so the table is not limited by the thing it replaces.

Needs mpmath, pinned to 1.3.0. The CI job that regenerates and diffs installs
exactly that version: this script exists to prove regeneration is
reproducible, and a library that changes its last digit would defeat that
quietly.
"""
import argparse

from mpmath import mp, mpf, quad, exp, sqrt, matrix, eigsy

X_MAX = mpf(55)          # above this the asymptotic form is exact in double
DEGREE = 13
WIDTH = mpf(5) / 2

# At n = 5 the measured switchover is NOT clean at X_MAX = 55: the largest
# root's weight, which decays fastest, is only within 1e-13 of the reference
# there (3.8e-13 at X = 56, falling to 2e-18 by X = 70 -- the asymptotic
# correction terms this form drops are themselves what take a few more units
# of X to become negligible). n <= 4 do not show this -- checked against the
# same reference, their X_MAX = 55 seam is already at the 1e-14 floor -- so
# only n = 5 moves its switchover out, to where the seam is clean again.
X_MAX_OVERRIDE = {5: mpf(65)}


def rys(n, X):
    """n Rys roots u_i and weights w_i for exp(-X u) on u in [0,1], u = t^2."""
    # moments M_k = int_0^1 t^(2k) exp(-X t^2) dt, i.e. the Boys function F_k(X)
    M = [quad(lambda t: t ** (2 * k) * exp(-X * t * t), [0, 1]) for k in range(2 * n + 1)]
    # Chebyshev algorithm: moments -> three-term recurrence coefficients
    sig_old = [mpf(0)] * (2 * n + 1)
    sig = list(M)
    a = [M[1] / M[0]]
    b = [mpf(0)]
    for k in range(1, n):
        new = [mpf(0)] * (2 * n + 1)
        for l in range(k, 2 * n - k + 1):
            new[l] = sig[l + 1] - a[k - 1] * sig[l] - (b[k - 1] * sig_old[l] if k > 1 else 0)
        a.append(new[k + 1] / new[k] - sig[k] / sig[k - 1])
        b.append(new[k] / sig[k - 1])
        sig_old, sig = sig, new
    # Golub-Welsch: eigenvalues of the Jacobi matrix are the roots, and the
    # weights are M_0 times the squared first components of the eigenvectors.
    J = matrix(n)
    for i in range(n):
        J[i, i] = a[i]
        if i:
            J[i, i - 1] = J[i - 1, i] = sqrt(b[i])
    u, V = eigsy(J)
    pairs = sorted([(u[i], M[0] * V[0, i] ** 2) for i in range(n)], key=lambda p: p[0])
    return [p[0] for p in pairs], [p[1] for p in pairs]


def cheb_nodes(lo, hi, deg):
    """The deg+1 Chebyshev nodes on [lo, hi], and the cosines the fit reuses."""
    N = deg + 1
    t = [mp.cos(mp.pi * (k + mpf(1) / 2) / N) for k in range(N)]
    xs = [mpf(1) / 2 * (hi - lo) * x + mpf(1) / 2 * (hi + lo) for x in t]
    cos_tab = [[mp.cos(mp.pi * j * (k + mpf(1) / 2) / N) for k in range(N)] for j in range(N)]
    return xs, cos_tab


def cheb_from_samples(ys, cos_tab):
    """Chebyshev coefficients from samples already taken at the nodes.

    Separate from the node generation so one `rys` evaluation per node serves
    all 2n series -- fitting each series independently would evaluate the same
    eigenproblem 2n times over.
    """
    N = len(ys)
    c = [2 * sum(ys[k] * cos_tab[j][k] for k in range(N)) / N for j in range(N)]
    c[0] /= 2
    return c


def fortran(x, digits=17):
    s = mp.nstr(x, digits, strip_zeros=False)
    if "e" not in s and "E" not in s and "." not in s:
        s += ".0"
    return s.replace("e", "e") + "_dp"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-n", "--nroots", type=int, required=True,
                    help="root count to tabulate: n covers every class with "
                         "L_total = la+lb+lc+ld in [2n-2, 2n-1]")
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("--dps", type=int, default=40)
    args = ap.parse_args()
    mp.dps = args.dps
    n = args.nroots
    x_max = X_MAX_OVERRIDE.get(n, X_MAX)

    nint = int(x_max / WIDTH)
    assert nint * WIDTH == x_max, "X_MAX must be a whole number of intervals"

    # coeff[interval][series][j]; series 0..n-1 are roots, n..2n-1 are weights
    coeff = []
    for iv in range(nint):
        lo, hi = iv * WIDTH, (iv + 1) * WIDTH
        xs, cos_tab = cheb_nodes(lo, hi, DEGREE)
        # one eigenproblem per node, shared by all 2n series
        samples = [rys(n, x if x > 0 else mpf("1e-30")) for x in xs]
        series = [cheb_from_samples([smp[0][i] for smp in samples], cos_tab) for i in range(n)]
        series += [cheb_from_samples([smp[1][i] for smp in samples], cos_tab) for i in range(n)]
        coeff.append(series)
        print(f"    n={n} interval {iv + 1}/{nint}  X in [{int(lo)}, {int(hi)}]")

    # Asymptotic constants, in closed form rather than by taking X huge.
    #
    # As X grows the measure exp(-X u) u^(-1/2) on [0,1] becomes, under
    # v = X u, the measure v^(-1/2) exp(-v) on [0, inf) -- Gauss-Laguerre with
    # alpha = -1/2, whose Jacobi matrix is known exactly:
    #     a_k = 2k + alpha + 1,   b_k = k (k + alpha)
    # Evaluating `rys` at a huge X instead loses the moments to underflow, the
    # recurrence coefficients go negative, and the eigensolver returns complex.
    alpha = mpf(-1) / 2
    J = matrix(n)
    for i in range(n):
        J[i, i] = 2 * i + alpha + 1
        if i:
            J[i, i - 1] = J[i - 1, i] = sqrt(mpf(i) * (i + alpha))
    cu, cV = eigsy(J)
    # total mass: the base rule carries u^(-1/2)/2, so the weights sum to
    # gamma(1/2)/2 = sqrt(pi)/2
    mass = sqrt(mp.pi) / 2
    pairs = sorted([(cu[i], mass * cV[0, i] ** 2) for i in range(n)], key=lambda q: q[0])
    a_r = [q[0] for q in pairs]
    a_w = [q[1] for q in pairs]

    flat = [c for iv in coeff for s in iv for c in s]
    lines = []
    for i in range(0, len(flat), 2):
        chunk = ", ".join(fortran(v) for v in flat[i:i + 2])
        lines.append("      " + chunk + ("," if i + 2 < len(flat) else ""))
    body = " &\n".join(lines)

    src = f'''!
! Rys roots and weights, {n} root{"s" if n != 1 else ""}, by piecewise Chebyshev interpolation.
!
! GENERATED by scripts/gen_rys_roots.py -n {n} -- do not edit.
!
! Degree {DEGREE} on intervals of width {float(WIDTH)} over X in [0, {int(x_max)}]; above that
! asymptotic form is exact in double precision. Reference values from mpmath at
! {args.dps} digits, so the table is not limited by the double precision
! construction (src/trc_rys.F90's docstring records a 1.3e-15 cross-check of
! that construction, which this table is not built from).
!
! LAYOUT: coefficient j of series s on interval iv is
!     RYS_C{n}(((iv-1)*{2*n} + (s-1))*{DEGREE+1} + j + 1)
! with series 1..{n} the roots and {n+1}..{2*n} the weights.
!
! NOT a parameter array: trc_boys records that nvfortran's fort2 segfaults on
! parameter arrays of this kind, and the same caution applies here.
!
module trc_rys_table{n}
   use trc_boys, only: dp
   implicit none
   private
   public :: RYS_NROOTS{n}, RYS_DEGREE{n}, RYS_NINT{n}, RYS_WIDTH{n}, RYS_XMAX{n}
   public :: RYS_C{n}, RYS_ASY_R{n}, RYS_ASY_W{n}, rys_table_init{n}

   integer,  parameter :: RYS_NROOTS{n} = {n}
   integer,  parameter :: RYS_DEGREE{n} = {DEGREE}
   integer,  parameter :: RYS_NINT{n}   = {nint}
   real(dp), parameter :: RYS_WIDTH{n}  = {float(WIDTH)}_dp
   real(dp), parameter :: RYS_XMAX{n}   = {int(x_max)}.0_dp

   real(dp), parameter :: RYS_ASY_R{n}({n}) = [ &
      {", ".join(fortran(v) for v in a_r)}]
   real(dp), parameter :: RYS_ASY_W{n}({n}) = [ &
      {", ".join(fortran(v) for v in a_w)}]

   real(dp) :: RYS_C{n}(RYS_NINT{n}*2*RYS_NROOTS{n}*(RYS_DEGREE{n} + 1)) = [ &
{body} &
      ]
   !$acc declare create(RYS_C{n})

contains

   !
   ! Push the table to the device. `declare create` allocates device storage
   ! but does NOT carry the host initialiser across -- a device copy full of
   ! garbage is exactly the silent-wrong-answer failure this project exists to
   ! avoid -- so the transfer is explicit and must be called once before any
   ! kernel. Same contract as boys_init.
   !
   subroutine rys_table_init{n}()
      !$acc update device(RYS_C{n})
   end subroutine rys_table_init{n}

end module trc_rys_table{n}
'''
    with open(args.output, "w") as fh:
        fh.write(src)
    print(f"  {args.output}: n={n}, {nint} intervals, {len(flat)} coefficients")


if __name__ == "__main__":
    main()
