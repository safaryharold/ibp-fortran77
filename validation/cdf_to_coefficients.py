"""Convert an official IBP coefficient CDF (SW_OPER_IBP_CLI_2__..._NNNN.cdf,
shipped inside the ibpmodel Python package) into the ASCII
ibp_coefficients.txt read by the FORTRAN 77 code.

Usage:  python cdf_to_coefficients.py <input.cdf> <output.txt>
Requires: numpy, cdflib.
"""
import sys
import numpy as np
import cdflib


def main(src, dst):
    with cdflib.CDF(src) as c:
        par = np.asarray(c.varget("Parameters"), float)        # C1 C2 TAU MU0 SIGMA
        g = np.asarray(c.varget("Intensity"), float)           # G0 GS1 GC1 GS2 GC2
        ts = np.asarray(c.varget("Monthly_LT_Shift"), float)   # (5, 12): T0 TS1 TC1 TS2 TC2
        phi = np.asarray(c.varget("Density_Estimators"), float)  # (360, 12)
        lons = np.asarray(c.varget("Density_Estimator_Lons"), float)
    assert par.shape == (5,) and g.shape == (5,) and ts.shape == (5, 12) and phi.shape == (360, 12)
    # The CDF labels the 360 longitude bins 1..360; the official code maps
    # bin k (0-based) to longitude -180 + k degrees, as the Fortran does.
    assert np.allclose(lons, np.arange(1, 361)), "unexpected longitude grid"
    row = lambda v: "".join(f"{x:25.16E}" for x in v)  # full double precision
    import os
    out = [f"IBP MODEL COEFFICIENTS (CHAMP/Swarm-derived), FORTRAN 77 data file, from {os.path.basename(src)}",
           "SCALAR PARAMS: C1 C2 TAU(=1/gamma,mean_lifetime_h) MU0 SIGMA", row(par),
           "DIURNAL OSC PARAMS: G0 GS(1) GC(1) GS(2) GC(2)", row(g),
           "MONTHLY LT-SHIFT T0MON (Jan..Dec, 12 values)", row(ts[0]),
           "MONTHLY HARMONIC-1 SIN COEFF TSMON(1,:) (Jan..Dec, 12 values)", row(ts[1]),
           "MONTHLY HARMONIC-1 COS COEFF TCMON(1,:) (Jan..Dec, 12 values)", row(ts[2]),
           "MONTHLY HARMONIC-2 SIN COEFF TSMON(2,:) (Jan..Dec, 12 values)", row(ts[3]),
           "MONTHLY HARMONIC-2 COS COEFF TCMON(2,:) (Jan..Dec, 12 values)", row(ts[4]),
           "LON-MONTH PHI TABLE: 360 rows (lon bin -180E..+179E, 1 deg step) x 12 cols (Jan..Dec)"]
    out += [row(r) for r in phi]
    with open(dst, "w") as f:
        f.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main(*sys.argv[1:3])
