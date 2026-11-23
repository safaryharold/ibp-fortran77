"""
plot_monthly_longitude_density.py

Reproduces the "Monthly Longitudinal Densities of EPD Occurrence
Probability" style figure: a 3x4 grid of panels (one per month),
each showing a weighted histogram of occurrence probability across
longitude (normalized to a density, i.e. summing to 1 rather than
showing raw probability), with a smoothed KDE curve overlaid.

"Weighted" here means: at each longitude bin, the model's peak-over-
local-time probability for that month is used as the WEIGHT of that
longitude bin, then the whole set of weights for the month is
normalized to a proper density. This shows *where* (which longitudes)
a given month's activity concentrates, not the absolute probability
level - that's a modeling choice made to match the reference figure's
apparent intent, since the F77 model itself has no built-in "density"
output.

Run with: python3 plot_monthly_longitude_density.py [loopresults.txt] [F10.7]
Requires: pandas, numpy, matplotlib, scipy
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

COLNAMES = ['DOY', 'Lon', 'F107', 'LT', 'Month', 'Prob']
MONTH_NAMES = ['JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN',
               'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']


def nearest(values, target):
    return min(values, key=lambda v: abs(v - target))


def weighted_hist_density(lon_values, weights, bins=25):
    """Weighted histogram normalized to a density (area = 1)."""
    counts, edges = np.histogram(lon_values, bins=bins, weights=weights,
                                  range=(-180, 180))
    widths = np.diff(edges)
    total = counts.sum() * widths.mean()
    density = counts / total if total > 0 else counts
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, density, edges, widths


def weighted_kde(lon_samples, weights, grid):
    """Gaussian KDE using longitude values repeated/weighted by
    resampling proportional to weight, since scipy's gaussian_kde
    takes a `weights` argument directly in recent versions."""
    kde = gaussian_kde(lon_samples, weights=weights)
    return kde(grid)


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'loopresults.txt'
    f107_target = float(sys.argv[2]) if len(sys.argv) > 2 else 120.0

    df = pd.read_csv(path, sep=r'\s+', skiprows=2, names=COLNAMES)
    f107_actual = min(df['F107'].unique(), key=lambda v: abs(v - f107_target))

    # Peak-over-LT probability per (DOY, Lon) at the chosen F10.7
    sub = df[np.isclose(df['F107'], f107_actual)]
    peak = sub.groupby(['DOY', 'Lon'])['Prob'].max().reset_index()
    peak = peak.merge(df[['DOY', 'Month']].drop_duplicates(), on='DOY')
    monthly = peak.groupby(['Month', 'Lon'])['Prob'].mean().reset_index()

    fig, axes = plt.subplots(3, 4, figsize=(20, 12))
    axes = axes.flatten()
    grid = np.linspace(-180, 180, 400)

    for ax, month in zip(axes, range(1, 13)):
        s = monthly[monthly['Month'] == month].sort_values('Lon')
        lon = s['Lon'].values
        w = s['Prob'].values

        centers, density, edges, widths = weighted_hist_density(lon, w,
                                                                   bins=25)
        ax.bar(centers, density, width=widths, color='#8fa8c8',
               edgecolor='none', alpha=0.9)

        if w.sum() > 0:
            kde_vals = weighted_kde(lon, w, grid)
            ax.plot(grid, kde_vals, color='red', linewidth=2)

        ax.set_title(MONTH_NAMES[month - 1], fontweight='bold')
        ax.set_xlim(-180, 180)
        if ax in (axes[0], axes[4], axes[8]):
            ax.set_ylabel('Density')
        if ax in axes[8:]:
            ax.set_xlabel('LON')

    fig.suptitle('Monthly Longitudinal Densities of IBP Occurrence '
                 f'Probability (F10.7 = {f107_actual:.0f} s.f.u)',
                 fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig('figures/Figure3_Monthly_Longitude_Density.png',
                dpi=150, bbox_inches='tight')
    print('Wrote figures/Figure3_Monthly_Longitude_Density.png')


if __name__ == '__main__':
    main()
