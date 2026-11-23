#!/usr/bin/env python3
"""
ibp_reference.py -- Independent Python reference implementation of the
Ionospheric Bubble Probability (IBP) model.

This is a self-contained re-implementation of the same published model the
FORTRAN 77 port implements (Stolle et al., 2024; official IBP release
v2.0.2, magnetic-field-based mode). It reads the identical plain-text coefficient file
(ibp_coefficients.txt) used by the Fortran code and evaluates the
occurrence probability with the same closed-form (error-function)
expression.

It serves two purposes in this supplement:
  1. An auditable, dependency-light cross-check of the Fortran output
     (see compare_504.py), independent of the official `ibpmodel`
     Python package.
  2. The model engine used by make_figures.py to regenerate the
     climatology figures.

Note on the error function: this reference uses Python's exact
`math.erf`, whereas the Fortran port uses the Abramowitz & Stegun 7.1.26
approximation `ERF77` (max abs error 1.5e-7). The small Fortran-minus-
reference residual reported by compare_504.py is therefore dominated by
that approximation, exactly as discussed in the paper.

The official model and coefficients live at
https://igit.iap-kborn.de/ibp/ibp-model (release v2.0.2; the coefficient
file here is converted from ..._0101.cdf by cdf_to_coefficients.py).
"""
import math
import os

_DEF_COEF = os.path.join(os.path.dirname(__file__), "ibp_coefficients.txt")


class IBPModel:
    def __init__(self, coeff_file=_DEF_COEF):
        self._load(coeff_file)

    def _load(self, path):
        with open(path) as fh:
            lines = [ln.rstrip("\n") for ln in fh]
        it = iter(lines)

        def row():
            # next non-header line -> list of floats; header lines are the
            # plain-text section titles, skipped exactly as READCF does.
            while True:
                ln = next(it)
                parts = ln.split()
                try:
                    return [float(p) for p in parts]
                except ValueError:
                    continue

        self.C1, self.C2, self.TAU, self.MU0, self.SIGMA = row()
        g = row()
        self.G0, self.GS1, self.GC1, self.GS2, self.GC2 = g
        self.T0 = row()
        self.TS1 = row()
        self.TC1 = row()
        self.TS2 = row()
        self.TC2 = row()
        self.PHI = [row() for _ in range(360)]
        self.LONST = -180.0
        self.NLON = 360

    @staticmethod
    def month_of_doy(doy):
        nd = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        t = doy
        for i, n in enumerate(nd, 1):
            if t <= n:
                return i
            t -= n
        return 0

    def g_osc(self, doy):
        r = doy / 365.0 * 2.0 * math.pi
        return (self.G0
                + self.GS1 * math.sin(r) + self.GC1 * math.cos(r)
                + self.GS2 * math.sin(2 * r) + self.GC2 * math.cos(2 * r))

    def phi_month(self, lon, month):
        a = math.fmod(lon - self.LONST, 360.0)
        if a < 0:
            a += 360.0
        idx = int(a / 1.0 + 1.0)
        idx = min(max(idx, 1), self.NLON)
        v = self.PHI[idx - 1][month - 1]
        return v if v > 0.0 else 0.0

    def ts_month(self, lon, month):
        r = lon / 360.0 * 2.0 * math.pi
        j = month - 1
        return (self.T0[j]
                + self.TS1[j] * math.sin(r) + self.TC1[j] * math.cos(r)
                + self.TS2[j] * math.sin(2 * r) + self.TC2[j] * math.cos(2 * r))

    def probability(self, doy, lon, f107, lt):
        """Occurrence probability P[X(t)=1], clamped to [0,1]."""
        month = self.month_of_doy(doy)
        if month < 1 or month > 12:          # invalid DOY -> default January
            month, doy = 1, 1
        t_int = math.fmod(lt + 12.0, 24.0)
        if t_int < 0:
            t_int += 24.0
        t_int -= 12.0

        lam1 = (self.g_osc(doy) + self.C1 + self.C2 * f107) * self.phi_month(lon, month)
        if lam1 < 0.0:
            lam1 = 0.0
        mu = self.MU0 + self.ts_month(lon, month)
        sig = self.SIGMA
        if self.TAU <= 1e-12 or sig <= 1e-12 or lam1 <= 1e-12:
            return 0.0
        lam2 = 1.0 / self.TAU            # survival decay rate = 1 / mean-lifetime

        omult = -1.0 / math.sqrt(2.0 * math.pi * sig * sig)
        xim = -1.0 / (2.0 * sig * sig)
        x0 = mu - lam2 / (2.0 * xim)
        y0 = lam2 * (mu - t_int - lam2 / (4.0 * xim))
        ims = math.sqrt(-xim)
        om2 = math.sqrt(math.pi) * omult * math.exp(y0) / ims
        xint = (0.5 + 0.5 * math.erf((t_int - x0) * ims)) * om2
        exparg = lam1 * xint
        if exparg < -700.0:
            prob = 1.0
        elif exparg > 1e-9:
            prob = 0.0
        else:
            prob = 1.0 - math.exp(exparg)
        return max(0.0, min(1.0, prob))


if __name__ == "__main__":
    import sys
    m = IBPModel()
    if len(sys.argv) == 5:
        doy, lon, f107, lt = (float(sys.argv[1]), float(sys.argv[2]),
                              float(sys.argv[3]), float(sys.argv[4]))
        print("P[X(t)=1] = %.7f" % m.probability(int(doy), lon, f107, lt))
    else:
        print("usage: python ibp_reference.py DOY LON F10.7 LT")
