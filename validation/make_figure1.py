"""Figure 1: data flow of the IBP model as implemented in the FORTRAN 77 port.

Writes Figure_1.pdf / Figure_1.eps (vector) and Figure_1.png (1000 dpi) to
the directory given as the first argument (default: figures/).
"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

plt.rcParams.update({"font.family": "DejaVu Sans", "mathtext.fontset": "dejavusans",
                     "font.size": 8.5, "pdf.fonttype": 42, "ps.fonttype": 42})

C = {  # muted, print- and colour-blind-safe fills; dark ink for text/edges
    "input": ("#e8eef5", "#4a6785"),
    "func":  ("#ffffff", "#4a6785"),
    "param": ("#fbefd9", "#a86b12"),
    "out":   ("#f6d88a", "#8a5a00"),
    "ink": "#1f2933", "coef": "#5b6670", "arrow": "#3d4852",
}

fig = plt.figure(figsize=(7.48, 4.3))          # 190 mm full ASR page width
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 58); ax.axis("off")

def box(x, y, w, h, title, sub=None, coef=None, kind="func", bold=False):
    fc, ec = C[kind]
    ax.add_patch(FancyBboxPatch((x - w/2, y - h/2), w, h, boxstyle="round,pad=0.25,rounding_size=1.2",
                                fc=fc, ec=ec, lw=1.0, zorder=2))
    lines = [(title, dict(fontsize=9, weight="bold" if bold else "normal", color=C["ink"]))]
    if sub:  lines.append((sub, dict(fontsize=7.3, color=C["ink"])))
    if coef: lines.append((coef, dict(fontsize=7.3, color=C["coef"], style="italic")))
    n = len(lines); step = 2.6
    for i, (s, kw) in enumerate(lines):
        ax.text(x, y + (n - 1) * step / 2 - i * step, s, ha="center", va="center", zorder=3, **kw)
    return dict(x=x, y=y, w=w, h=h)

def side(b, s, d=0.0):
    pad = 0.35
    return {"l": (b["x"] - b["w"]/2 - pad, b["y"] + d), "r": (b["x"] + b["w"]/2 + pad, b["y"] + d),
            "t": (b["x"] + d, b["y"] + b["h"]/2 + pad), "b": (b["x"] + d, b["y"] - b["h"]/2 - pad)}[s]

def route(*pts):
    """Orthogonal polyline through pts; arrowhead on the last segment."""
    xs, ys = zip(*pts)
    if len(pts) > 2:
        ax.plot(xs[:-1], ys[:-1], color=C["arrow"], lw=0.9, zorder=1, solid_capstyle="butt")
    ax.add_patch(FancyArrowPatch(pts[-2], pts[-1], arrowstyle="-|>", mutation_scale=9, lw=0.9,
                                 color=C["arrow"], shrinkA=0, shrinkB=0, zorder=1))

def hv(p, q, xm=None):
    """Horizontal-vertical-horizontal elbow from p to q, turning at x = xm."""
    if abs(p[1] - q[1]) < 1e-9: return route(p, q)
    xm = (p[0] + q[0]) / 2 if xm is None else xm
    route(p, (xm, p[1]), (xm, q[1]), q)

def lane(x0, x1, label):
    ax.text((x0 + x1)/2, 55.6, label, ha="center", va="center", fontsize=8.2, weight="bold", color=C["coef"])
    ax.plot([x0 + 1, x1 - 1], [54.0, 54.0], color="#c3cad2", lw=0.8)

lane(0, 21, "Inputs"); lane(21, 50, "Model functions"); lane(50, 74, "Physical parameters"); lane(74, 100, "Probability")

# inputs
f107 = box(10.5, 49,   17, 5, "Solar flux", "F10.7 (sfu)", kind="input")
doy  = box(10.5, 40,   17, 5, "Day of year", "doy (1–365)", kind="input")
mon  = box(10.5, 31,   17, 4.2, "month", "derived from doy", kind="input")
lon  = box(10.5, 21,   17, 5, "Longitude", "lon (°E)", kind="input")
lt   = box(10.5, 7,    17, 5, "Local time", "LT (h)", kind="input")

# model functions
gosc = box(35.5, 44,   24, 7.2, r"$g_{\mathrm{osc}}(\mathrm{doy})$", "seasonal harmonics, Eq. (4)",
           r"$g_0,\ g_i^{(s)},\ g_i^{(c)}$")
phi  = box(35.5, 29.5, 24, 7.2, r"$\varphi_{\mathrm{month}}(\mathrm{lon})$", "kernel-density lookup table",
           "1° longitude × 12 months")
ts   = box(35.5, 17,   24, 7.2, r"$\mathrm{ts}_{\mathrm{month}}(\mathrm{lon})$", "time-shift harmonics, Eq. (6)",
           r"$t_0,\ t_i^{(s)},\ t_i^{(c)}$ per month")
tw   = box(35.5, 7,    24, 5.4, r"$t$ (h from local midnight)", "LT wrapped to [−12, 12)", None)

# parameters
lam  = box(62, 38,   21, 9.2, r"$\lambda$", "bubble intensity, Eq. (3)",
           r"$[g_{\mathrm{osc}}+C_1+C_2\,\mathrm{F10.7}]\,\varphi_{\mathrm{month}}$", kind="param", bold=True)
glob = box(62, 27.2, 21, 4.6, r"$\sigma,\ \ \tau = 1/\gamma$", "global constants", None, kind="param")
mu   = box(62, 17,   21, 7.2, r"$\mu$", "mean onset time, Eq. (5)", r"$\mu_0+\mathrm{ts}_{\mathrm{month}}$",
           kind="param", bold=True)

# probability
I    = box(87, 17,   22, 9.2, r"$I(t,\gamma,\mu,\sigma)$", "onset–survival integral, Eq. (2)",
           "closed form with erf", kind="func")
P    = box(87, 38,   22, 9.2, r"$P = 1-\exp(\lambda\, I)$", "occurrence probability, Eq. (1)", r"$0 \leq P \leq 1$",
           kind="out", bold=True)

def dot(x, y): ax.plot([x], [y], "o", ms=2.6, color=C["arrow"], zorder=1)

# intensity path
hv(side(doy, "r"), side(gosc, "l"), xm=21.5)
route(side(doy, "b"), side(mon, "t"))
p = side(f107, "r"); q = side(lam, "t"); route(p, (q[0], p[1]), q)
hv(side(gosc, "r"), side(lam, "l", 2.4), xm=49.5)
hv(side(phi, "r"), side(lam, "l", -2.4), xm=49.5)
# (month, lon) bus feeding both monthly functions
xb = 21.5
ax.plot([side(mon, "r")[0], xb], [mon["y"]] * 2, color=C["arrow"], lw=0.9, zorder=1)
ax.plot([side(lon, "r")[0], xb], [lon["y"]] * 2, color=C["arrow"], lw=0.9, zorder=1)
ax.plot([xb, xb], [mon["y"], ts["y"]], color=C["arrow"], lw=0.9, zorder=1)
route((xb, phi["y"]), side(phi, "l")); route((xb, ts["y"]), side(ts, "l"))
for y in (mon["y"], lon["y"], phi["y"]): dot(xb, y)
# timing path
route(side(ts, "r"), side(mu, "l"))
route(side(lt, "r"), side(tw, "l"))
route(side(mu, "r"), side(I, "l"))
p = side(tw, "r"); q = side(I, "l", -2.8); route(p, (74.5, p[1]), (74.5, q[1]), q)
p = side(glob, "r"); q = side(I, "l", 2.8); route(p, (74.5, p[1]), (74.5, q[1]), q)
# combine
route(side(lam, "r"), side(P, "l"))
route(side(I, "t"), side(P, "b"))

ax.text(50, 0.6, r"Grey italic lines give the coefficients or formula used. All coefficients, including the "
        r"$\varphi_{\mathrm{month}}$ table, are read from ibp_coefficients.txt (Tables 1–2).",
        ha="center", va="center", fontsize=6.8, color=C["coef"])

out = sys.argv[1] if len(sys.argv) > 1 else "figures"
os.makedirs(out, exist_ok=True)
for ext, kw in (("pdf", {}), ("eps", {}), ("png", {"dpi": 1000})):
    fig.savefig(os.path.join(out, f"Figure_1.{ext}"), **kw)
print("written to", out)
