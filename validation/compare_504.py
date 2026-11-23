#!/usr/bin/env python3
"""
compare_504.py -- Reproduce the Fortran-vs-reference cross-validation.

Reads the 504-case input grid (batch_inputs_504.txt), evaluates the
occurrence probability with the independent Python reference
(ibp_reference.py), and compares it row-by-row against the committed
Fortran output (fortran_loopresults_504.txt). Reports the maximum and
mean absolute probability difference.

The residual is dominated by the Fortran port's ERF77 error-function
approximation (the reference uses exact math.erf); see ibp_reference.py.

Run:
    python compare_504.py
"""
import os
from ibp_reference import IBPModel

HERE = os.path.dirname(os.path.abspath(__file__))


def read_inputs(path):
    rows = []
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 4:
                continue
            try:
                rows.append((int(float(p[0])), float(p[1]), float(p[2]), float(p[3])))
            except ValueError:
                continue            # header line
    return rows


def read_fortran(path):
    probs = []
    with open(path) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) < 6:
                continue
            try:
                int(p[0])
            except ValueError:
                continue            # header / rule line
            probs.append(float(p[-1]))
    return probs


def main():
    model = IBPModel(os.path.join(HERE, "ibp_coefficients.txt"))
    inputs = read_inputs(os.path.join(HERE, "batch_inputs_504.txt"))
    fort = read_fortran(os.path.join(HERE, "fortran_loopresults_504.txt"))
    if len(inputs) != len(fort):
        raise SystemExit("row count mismatch: %d inputs vs %d fortran rows"
                         % (len(inputs), len(fort)))

    max_abs = 0.0
    sum_abs = 0.0
    worst = None
    for (doy, lon, f107, lt), pf in zip(inputs, fort):
        pr = model.probability(doy, lon, f107, lt)
        d = abs(pr - pf)
        sum_abs += d
        if d > max_abs:
            max_abs, worst = d, (doy, lon, f107, lt, pr, pf)
    n = len(fort)
    print("cases compared        : %d" % n)
    print("max |P_ref - P_fortran|: %.3e" % max_abs)
    print("mean |P_ref - P_fortran|: %.3e" % (sum_abs / n))
    if worst:
        print("worst case            : DOY=%d lon=%.1f F10.7=%.1f LT=%.2f  "
              "ref=%.7f fortran=%.7f" % worst)


if __name__ == "__main__":
    main()
