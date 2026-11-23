"""Figure 8: component-level and full-domain comparison of the FORTRAN 77
code with the official IBP Python package.

(a) For every intermediate model quantity -- g_osc(doy), ts_month(lon),
    phi_month(lon), lambda, mu -- and for the error function and the final
    probability P, the maximum Fortran-minus-official difference over the
    same 50,000 random inputs used by compare_official.py, normalised by the
    largest magnitude of that quantity.
(b) Fortran P against official P for those inputs, with the distribution of
    the differences as an inset.

Requirements: numpy, scipy, matplotlib, cdflib and the official package
    pip install ibpmodel==2.0.2
(ibpmodel 2.0 imports PyTorch for its optional neural-network mode, which
is not used here). Build the component driver first (see ibp_components.f),
then run
    python make_figure8.py
Writes figures/figure8_components.png and .pdf and prints the numbers.
"""
import os
import subprocess

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy import special
from ibpmodel import ibpcalc

try:
    import seaborn as sns
    sns.set_theme(context="notebook", style="whitegrid")
except ImportError:
    pass
mpl.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 300,
    "font.family": "serif", "mathtext.fontset": "cm",
    "font.size": 14, "axes.titlesize": 16, "axes.titleweight": "bold",
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 12, "axes.linewidth": 1.0,
    "savefig.bbox": "tight", "pdf.fonttype": 42, "savefig.pad_inches": 0.15,
})

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
EPS = np.finfo(float).eps


def main():
    # same random inputs as compare_official.py
    rng = np.random.default_rng(7)
    n = 50000
    doy = rng.integers(1, 366, n)
    lon = rng.integers(-180, 180, n).astype(float)
    f107 = rng.uniform(80, 200, n)
    lt = rng.uniform(-6, 24, n)
    xerf = np.linspace(-6.0, 6.0, n)            # arguments for the erf check

    inp = "\n".join(f"{d} {l:.1f} {f:.10f} {t:.10f} {x:.15e}"
                    for d, l, f, t, x in zip(doy, lon, f107, lt, xerf))
    out = subprocess.run([os.path.join(HERE, "ibp_components")], input=inp,
                         capture_output=True, text=True, cwd=HERE, check=True).stdout
    F = np.array([[float(v) for v in r.split()] for r in out.splitlines()])
    f_mon, f_g, f_ts, f_phi, f_lam, f_mu, f_p, f_erf = F.T
    # use exactly the values the Fortran program read
    f107 = np.array([float(f"{v:.10f}") for v in f107])
    lt = np.array([float(f"{v:.10f}") for v in lt])
    xerf = np.array([float(f"{v:.15e}") for v in xerf])

    # official ibpmodel quantities
    data = ibpcalc.read_model_file()
    _, month = ibpcalc.align_time_of_year(doy, np.zeros_like(doy))
    params = np.asarray(data["Parameters"])
    dens = np.asarray(data["Density_Estimators"])
    o_g = ibpcalc.fourier_model(np.asarray(data["Intensity"]), doy, 365.0)
    o_ts = ibpcalc.fourier_model(np.asarray(data["Monthly_LT_Shift"])[:, month - 1], lon, 360.0)
    o_phi = dens.reshape(-1)[np.array((lon + 180.0) * 12 + month - 1, dtype=int)]
    o_lam = ibpcalc.compute_lambda(lon, params, f107, o_g, dens, month - 1)
    o_mu = params[3] + o_ts
    o_p = ibpcalc.compute_probability_exp(doy, np.zeros_like(doy), lon, lt, f107, data)
    o_erf = special.erf(xerf)

    assert np.array_equal(f_mon.astype(int), month), "month mapping differs"

    rows = [(r"$g_{\mathrm{osc}}$", f_g, o_g), (r"$\mathrm{ts}_{\mathrm{month}}$", f_ts, o_ts),
            (r"$\varphi_{\mathrm{month}}$", f_phi, o_phi), (r"$\lambda$", f_lam, o_lam),
            (r"$\mu$", f_mu, o_mu), ("erf (ERF77)", f_erf, o_erf), (r"$P$", f_p, o_p)]
    labels, vals = [], []
    print("quantity            max|diff|        max|x|        normalised")
    for name, a, b in rows:
        dmax = np.abs(a - b).max(); scale = np.abs(b).max()
        labels.append(name); vals.append(dmax / scale)
        print(f"{name:20s} {dmax:12.3e} {scale:12.3e} {dmax/scale:12.3e}")
    dp = f_p - o_p
    print(f"P: max|diff| {np.abs(dp).max():.2e}  mean|diff| {np.abs(dp).mean():.2e}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6.2), constrained_layout=True,
                                   gridspec_kw={"width_ratios": [1.05, 1]})
    # ---- (a) component-wise differences
    floor = 1e-18
    colors = ["#4a6785"] * 5 + ["#a86b12", "#8a5a00"]
    y = np.arange(len(labels))[::-1]
    ax1.barh(y, [v if v > 0 else 0 for v in vals], color=colors, height=0.6, left=floor)
    ax1.set_xscale("log"); ax1.set_xlim(floor, 1e-4)
    ax1.set_yticks(y); ax1.set_yticklabels(labels)
    ax1.axvline(EPS, color="0.35", ls="--", lw=1.2,
                label=r"double-precision rounding ($2.2\times10^{-16}$)")
    ax1.legend(loc="upper right", fontsize=11, frameon=True, framealpha=0.95)
    def sci(v):
        m, e = f"{v:.1e}".split("e")
        return rf"${m}\times10^{{{int(e)}}}$"
    for yi, v in zip(y, vals):
        if v == 0:
            ax1.text(floor * 3, yi, "identical (difference = 0)", va="center", fontsize=11)
        else:
            ax1.text(v * 1.8, yi, sci(v), va="center", fontsize=11)
    ax1.set_xlabel("Max |Fortran − official| / max |official|")
    ax1.set_title("(a)", loc="left")
    ax1.grid(axis="y", visible=False)

    # ---- (b) Fortran vs official probability, all inputs
    hb = ax2.hexbin(o_p, f_p, gridsize=70, cmap="rocket_r" if "rocket_r" in plt.colormaps() else "magma_r",
                    norm=LogNorm(), mincnt=1, extent=(0, 1, 0, 1))
    ax2.plot([0, 1], [0, 1], color="0.2", lw=0.9, ls="--")
    ax2.set_xlim(0, 1); ax2.set_ylim(0, 1); ax2.set_aspect("equal")
    ax2.set_xlabel("Official ibpmodel 2.0.2 probability")
    ax2.set_ylabel("FORTRAN 77 probability")
    ax2.set_title("(b)", loc="left")
    cb = fig.colorbar(hb, ax=ax2, pad=0.02, shrink=0.9); cb.set_label("Number of inputs")
    ins = ax2.inset_axes([0.50, 0.08, 0.46, 0.34])
    ins.hist(dp * 1e6, bins=80, color="#4a6785")
    ins.set_yscale("log"); ins.tick_params(labelsize=9)
    ins.set_xlabel(r"Fortran − official ($\times 10^{-6}$)", fontsize=9.5)
    ins.set_ylabel("Count", fontsize=9.5)
    ins.set_facecolor("white")
    for s in ins.spines.values(): s.set_edgecolor("0.4")

    os.makedirs(OUT, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"figure8_components.{ext}"))
    print("wrote figures/figure8_components.png/.pdf")


if __name__ == "__main__":
    main()
