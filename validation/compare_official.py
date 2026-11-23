"""Compare the FORTRAN 77 library against the official IBP Python package.

Draws 50,000 random inputs (fixed seed), evaluates them with the Fortran
library through ibp_dump and with the official ibpmodel code, and reports
the differences with the exact erf and with the ERF77 approximation used
by the Fortran code substituted into the official code.

Requirements: numpy, scipy, cdflib and the official package
    pip install ibpmodel==2.0.2   (the coefficients shipped with this code)
(ibpmodel 2.0 also imports PyTorch for its optional neural-network mode,
which is not used here.)
Build ibp_dump first (see ibp_dump.f), then run
    python compare_official.py [coefficient.cdf]
Without an argument the CDF bundled with ibpmodel is used. To test another
release, convert its CDF with cdf_to_coefficients.py into
ibp_coefficients.txt next to ibp_dump, and pass the same CDF here.
"""
import os, subprocess, sys
import numpy as np
from ibpmodel import ibpcalc

HERE = os.path.dirname(os.path.abspath(__file__))


def erf77(x):  # Abramowitz & Stegun 7.1.26, as in pbutils.f
    ax = np.abs(x); t = 1.0 / (1.0 + 0.3275911 * ax)
    y = 1.0 - (((((1.061405429 * t - 1.453152027) * t) + 1.421413741) * t
                - 0.284496736) * t + 0.254829592) * t * np.exp(-ax * ax)
    return np.sign(x) * y


def main():
    cdf = sys.argv[1] if len(sys.argv) > 1 else None
    data = ibpcalc.read_model_file(cdf)
    rng = np.random.default_rng(7)
    n = 50000
    doy = rng.integers(1, 366, n)
    lon = rng.integers(-180, 180, n).astype(float)
    f107 = rng.uniform(80, 200, n)
    lt = rng.uniform(-6, 24, n)
    inp = "\n".join(f"{d} {l:.1f} {f:.10f} {t:.10f}" for d, l, f, t in zip(doy, lon, f107, lt))
    out = subprocess.run([os.path.join(HERE, "ibp_dump")], input=inp, capture_output=True,
                         text=True, cwd=HERE, check=True).stdout
    pf = np.array([float(r.split()[4]) for r in out.splitlines()])
    special = ibpcalc.special
    exact = special.erf
    try:
        p_exact = ibpcalc.compute_probability_exp(doy, np.zeros_like(doy), lon, lt, f107, data)
        special.erf = erf77
        p_erf77 = ibpcalc.compute_probability_exp(doy, np.zeros_like(doy), lon, lt, f107, data)
    finally:
        special.erf = exact
    d1 = np.abs(pf - p_exact); d2 = np.abs(pf - p_erf77)
    print(f"cases                                 : {n}")
    print(f"Fortran vs official (exact erf)       : max {d1.max():.2e}  mean {d1.mean():.2e}")
    print(f"Fortran vs official (ERF77 substituted): max {d2.max():.2e}")


if __name__ == "__main__":
    main()
