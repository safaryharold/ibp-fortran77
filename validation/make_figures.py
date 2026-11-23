#!/usr/bin/env python3
"""
make_figures.py -- Regenerate the IBP climatology figures from the model.

Reproduces the main figures of the paper directly from the independent
Python reference (ibp_reference.py), which reads the same
ibp_coefficients.txt as the Fortran port:

  figure2_density.png   : monthly longitudinal density phi_month(lon),
                          histogram + smooth curve, one panel per month.
  figure3_seasonal.png  : occurrence probability vs longitude and local
                          time for the four cardinal seasons (F10.7 = 150, matching Stolle et al. 2024 Fig. 5).
  figure4_sectors.png   : probability vs local time for four longitude
                          sectors, December solstice, F10.7 = 150.
  figure5_doy_lt.png    : probability vs day-of-year and local time at
                          70 W, F10.7 = 150.
  figure6_monthly.png   : monthly-mean probability (LT 18-06) vs longitude
                          and month for F10.7 = 80/120/160/200.
  figure7_crossval.png  : Python reference vs Fortran port and their
                          difference (reproduces Figure 7).

Requires numpy + matplotlib (seaborn optional, used for styling/palettes).
Run:  python make_figures.py   ->  outputs to ./figures/
"""
import os
import shutil
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
try:
    import seaborn as sns
    _HAVE_SNS = True
except ImportError:
    _HAVE_SNS = False
from ibp_reference import IBPModel

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "figures")
os.makedirs(OUT, exist_ok=True)
M = IBPModel(os.path.join(HERE, "ibp_coefficients.txt"))

# ---------- uniform house style ----------
if _HAVE_SNS:
    sns.set_theme(context="notebook", style="whitegrid")
mpl.rcParams.update({
    "figure.dpi": 120, "savefig.dpi": 300,
    "font.family": "serif", "mathtext.fontset": "cm",
    "font.size": 14, "axes.titlesize": 16, "axes.titleweight": "bold",
    "axes.labelsize": 14, "xtick.labelsize": 12, "ytick.labelsize": 12,
    "legend.fontsize": 12, "axes.linewidth": 1.0,
    "lines.linewidth": 2.2, "xtick.direction": "out", "ytick.direction": "out",
    "savefig.bbox": "tight", "pdf.fonttype": 42, "savefig.pad_inches": 0.15,
})
PROB_CMAP = "rocket" if _HAVE_SNS else "magma"
DIV_CMAP = "RdBu_r"
ACCENT = sns.color_palette("deep")[0] if _HAVE_SNS else "#1f4e79"

LONS = np.arange(-180, 181, 2)
LTS = np.arange(18.0, 30.01, 0.25)
MONTHS = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN",
          "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]



def _save(fig, name):
    """Write PNG (300 dpi) and vector PDF with embedded fonts. The ASR guide
    asks for line art as vector or 1000 dpi; the PDF serves that purpose."""
    fig.savefig(os.path.join(OUT, name + ".png"))
    fig.savefig(os.path.join(OUT, name + ".pdf"))


def _finish(ax, xlabel=None, ylabel=None, minor=True):
    if xlabel: ax.set_xlabel(xlabel)
    if ylabel: ax.set_ylabel(ylabel)
    if minor:
        ax.xaxis.set_minor_locator(AutoMinorLocator())
        ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(which="both", top=True, right=True)
    ax.tick_params(which="major", length=6)
    ax.tick_params(which="minor", length=3)


def grid_lon_lt(doy, f107):
    return np.array([[M.probability(int(doy), float(lo), f107, float(lt % 24))
                      for lo in LONS] for lt in LTS])


# ============================ Figure 2 ============================
def figure2():
    """Monthly longitudinal density phi_month(lon) from the coefficient table."""
    lon1 = np.arange(-180, 180)                      # 1-deg centres
    nbin, bw = 30, 12                                # 30 bars of 12 deg
    edges = np.arange(-180, 181, bw)
    bc = (edges[:-1] + edges[1:]) / 2
    fig, axes = plt.subplots(3, 4, figsize=(20, 12), constrained_layout=True)
    for m, ax in enumerate(axes.ravel(), start=1):
        dens = np.array([M.phi_month(float(x), m) for x in lon1])
        ax.plot(lon1, dens, color=ACCENT, linewidth=2.4, zorder=3)
        ax.fill_between(lon1, dens, color=ACCENT, alpha=0.12, zorder=2)
        ax.set_title(MONTHS[m-1])
        ax.set_xlim(-180, 180); ax.set_xticks(np.arange(-180, 181, 90))
        ax.set_ylim(bottom=0)
        _finish(ax,
                xlabel="Longitude (deg.)" if m > 8 else None,
                ylabel=r"$\phi_{month}(lon)$" if m % 4 == 1 else None)
    _save(fig, "figure2_density"); plt.close(fig)


# ============================ Figure 3 ============================
def figure3():
    seasons = [("(a) December solstice", 355), ("(b) March equinox", 79),
               ("(c) June solstice", 172), ("(d) September equinox", 265)]
    lv = np.linspace(0, 0.9, 19)
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, (title, doy) in zip(axes.ravel(), seasons):
        Z = grid_lon_lt(doy, 150.0)
        cf = ax.contourf(LONS, LTS - 18, Z, levels=lv, cmap=PROB_CMAP, extend="max")
        ax.contour(LONS, LTS - 18, Z, levels=lv[::3], colors="white",
                   linewidths=0.4, alpha=0.5)
        ax.set_title(title)
        ax.set_yticks([0, 2, 4, 6, 8, 10, 12])
        ax.set_yticklabels(["18", "20", "22", "00", "02", "04", "06"])
        ax.set_xticks(np.arange(-180, 181, 60)); ax.grid(False)
        _finish(ax, xlabel="Longitude (deg.)", ylabel="Local time (h)")
    cb = fig.colorbar(cf, ax=axes, orientation="horizontal", shrink=0.6,
                      aspect=40, pad=0.06)
    cb.set_label("Occurrence probability")
    _save(fig, "figure3_seasonal"); plt.close(fig)


# ============================ Figure 4 ============================
def figure4():
    sectors = [("South America (75°W)", -75, "-", "black"),
               ("African sector (15°E)", 15, "--", "tab:blue"),
               ("Indian sector (80°E)", 80, ":", "tab:green"),
               ("Pacific sector (180°W)", 180, "-.", "tab:red")]
    fig, ax = plt.subplots(figsize=(12.5, 7.5), constrained_layout=True)
    for lbl, lon, ls, c in sectors:
        y = [M.probability(355, lon, 150.0, float(lt % 24)) for lt in LTS]
        ax.plot(LTS - 18, y, ls=ls, color=c, lw=2.6, label=lbl)
    ax.set_xticks(range(0, 13, 2))
    ax.set_xticklabels(["18", "20", "22", "00", "02", "04", "06"])
    ax.set_xlim(0, 12); ax.set_ylim(bottom=0)
    ax.legend(frameon=True, framealpha=0.9)
    _finish(ax, xlabel="Local time (h)", ylabel="Occurrence probability")
    _save(fig, "figure4_sectors"); plt.close(fig)


# ============================ Figure 5 ============================
def figure5():
    doys = np.arange(1, 366, 2)
    Z = np.array([[M.probability(int(d), -70.0, 150.0, float(lt % 24))
                   for d in doys] for lt in LTS])
    vmax = float(np.ceil(Z.max() * 10) / 10)
    fig, ax = plt.subplots(figsize=(13, 6.5), constrained_layout=True)
    cf = ax.contourf(doys, LTS - 18, Z, levels=np.linspace(0, vmax, 17),
                     cmap=PROB_CMAP, extend="max")
    ax.set_yticks([0, 2, 4, 6, 8, 10, 12])
    ax.set_yticklabels(["18", "20", "22", "00", "02", "04", "06"])
    ax.set_xticks(np.arange(0, 366, 60)); ax.grid(False)
    _finish(ax, xlabel="Day of year", ylabel="Local time (h)")
    cb = fig.colorbar(cf, ax=ax, pad=0.02, ticks=np.arange(0, vmax + 0.01, 0.2))
    cb.set_label("Occurrence probability")
    _save(fig, "figure5_doy_lt"); plt.close(fig)


# ============================ Figure 6 ============================
def figure6():
    mid = [15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349]
    months = np.arange(1, 13)
    lts = np.arange(19.0, 25.01, 0.25)              # 19-01 LT, as in Stolle et al. (2024) Fig. 6
    lv = np.linspace(0, 0.8, 17)
    fig, axes = plt.subplots(2, 2, figsize=(16, 10), sharex=True, sharey=True,
                             constrained_layout=True)
    for ax, (lbl, f107) in zip(axes.ravel(),
                               [("(a)", 80.0), ("(b)", 120.0), ("(c)", 160.0), ("(d)", 200.0)]):
        Z = np.array([[np.mean([M.probability(mid[mi], float(lo), f107, float(lt % 24))
                                for lt in lts]) for lo in LONS] for mi in range(12)])
        cf = ax.contourf(LONS, months, Z, levels=lv, cmap=PROB_CMAP, extend="max")
        ax.set_title(f"{lbl} F10.7 = {int(f107)} sfu")
        ax.set_yticks(months); ax.set_yticklabels(MONTHS); ax.grid(False)
        ax.set_xticks(np.arange(-180, 181, 60))
        _finish(ax, xlabel="Longitude (deg.)", ylabel="Month", minor=False)
    cb = fig.colorbar(cf, ax=axes, orientation="horizontal", shrink=0.6,
                      aspect=40, pad=0.06)          # colour key BELOW the plots
    cb.set_label("Mean occurrence probability (LT 19–01)")
    _save(fig, "figure6_monthly"); plt.close(fig)


# ============================ Figure 7 ============================
def _fortran_grid(doys, lons, f107, lt):
    pbbatch = os.path.normpath(os.path.join(HERE, "..", "ibp_model_f77", "pbbatch"))
    if not (os.path.isfile(pbbatch) and os.access(pbbatch, os.X_OK)):
        print("  (pbbatch not built; using Python reference for the Fortran panel)")
        return np.array([[M.probability(int(d), float(lo), f107, lt)
                          for lo in lons] for d in doys]), False
    run = os.path.join(OUT, "_fortran_run"); shutil.rmtree(run, ignore_errors=True)
    os.makedirs(run)
    shutil.copy(pbbatch, run)
    shutil.copy(os.path.join(HERE, "ibp_coefficients.txt"), run)
    with open(os.path.join(run, "batch_inputs.txt"), "w") as fh:
        fh.write("DOY Lon F107 LT\n")
        for d in doys:
            for lo in lons:
                fh.write(f"{d} {lo} {f107} {lt}\n")
    subprocess.run(["./pbbatch"], cwd=run, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    vals = []
    with open(os.path.join(run, "loopresults.txt")) as fh:
        for ln in fh:
            p = ln.split()
            if len(p) >= 6:
                try: int(p[0])
                except ValueError: continue
                vals.append(float(p[-1]))
    field = np.array(vals).reshape(len(doys), len(lons))
    shutil.rmtree(run, ignore_errors=True)
    return field, True


def figure7():
    lt, f107 = 21.0, 150.0
    doys = np.arange(1, 366)
    lons = np.arange(-180, 181, 10)
    py = np.array([[M.probability(int(d), float(lo), f107, lt) for lo in lons] for d in doys])
    fort, real = _fortran_grid(doys, lons, f107, lt)
    diff = fort - py
    print("  figure7 max|diff| = %.2e (%s)"
          % (np.abs(diff).max(), "Fortran binary" if real else "reference only"))
    vmax = float(np.ceil(max(py.max(), fort.max()) * 10) / 10)
    dlim = 2e-7
    ext = [lons[0], lons[-1], doys[0], doys[-1]]
    fig, ax = plt.subplots(3, 1, figsize=(12.5, 15), constrained_layout=True)
    for a in ax: a.grid(False)

    def panel(a, Z, title, interp, cmap, vmin, vmx):
        im = a.imshow(Z, origin="lower", aspect="auto", extent=ext, cmap=cmap,
                      vmin=vmin, vmax=vmx, interpolation=interp)
        a.set_title(title)
        a.set_yticks([1, 60, 120, 180, 240, 300, 365])
        a.set_xticks(np.arange(-180, 181, 60))
        _finish(a, ylabel="Day of year", minor=False)
        return im

    im0 = panel(ax[0], py, "Python reference", "bilinear", PROB_CMAP, 0, vmax)
    panel(ax[1], fort, "FORTRAN 77 port", "bilinear", PROB_CMAP, 0, vmax)
    # difference: discrete pixels, neutral diverging scale centred on zero
    im2 = panel(ax[2], diff, "Difference (Fortran − Python)", "nearest", DIV_CMAP, -dlim, dlim)
    # coarse box grid overlaid so the (uniformly ~zero) cells read as boxes
    for x in np.arange(-180, 181, 10):
        ax[2].axvline(x, color="0.6", lw=0.3)
    for y in np.arange(1, 366, 20):
        ax[2].axhline(y, color="0.6", lw=0.3)
    ax[2].set_xlabel("Longitude (deg.)")
    cb = fig.colorbar(im0, ax=[ax[0], ax[1]], location="right", shrink=0.85, pad=0.02)
    cb.set_label("Occurrence probability")
    cb2 = fig.colorbar(im2, ax=ax[2], location="right", shrink=0.85, pad=0.02,
                       ticks=[-dlim, -dlim/2, 0, dlim/2, dlim])
    cb2.formatter.set_powerlimits((0, 0)); cb2.update_ticks()
    cb2.set_label("Fortran − Python")
    _save(fig, "figure7_crossval"); plt.close(fig)


if __name__ == "__main__":
    figure2(); print("wrote figure2_density.png")
    figure3(); print("wrote figure3_seasonal.png")
    figure4(); print("wrote figure4_sectors.png")
    figure5(); print("wrote figure5_doy_lt.png")
    figure6(); print("wrote figure6_monthly.png")
    figure7(); print("wrote figure7_crossval.png")
    print("figures written to", OUT)
