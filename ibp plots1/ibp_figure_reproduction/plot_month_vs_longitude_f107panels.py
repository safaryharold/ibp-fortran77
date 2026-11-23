"""
plot_month_vs_longitude_f107panels.py

Reproduces the "IBP Model: Occurrence Probability vs. Month and
Longitude" style figure: a 2x2 grid of filled-contour panels, one per
F10.7 level, with Longitude on the x-axis and Month on the y-axis
(Jan at the bottom, Dec at the top). Each panel's top edge is labeled
with rough longitude-sector names, matching the reference figure.

The plotted quantity is the PEAK occurrence probability over local
time for each (month, longitude) pair - the model has no single
"index" of its own, so peak-over-LT is used as a representative
per-night summary (documented here since it's a modeling choice, not
something the F77 program outputs directly).

Run with: python3 plot_month_vs_longitude_f107panels.py [loopresults.txt]
Requires: pandas, numpy, matplotlib
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

COLNAMES = ['DOY', 'Lon', 'F107', 'LT', 'Month', 'Prob']
MONTH_NAMES = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# Approximate longitude-sector labels for the secondary top axis,
# matching common ionospheric-literature naming. Edit freely.
SECTOR_LABELS = {-180: 'Pacific', -90: 'America', 0: 'Atlantic Africa',
                  80: 'India', 180: 'Pacific'}

F107_PANELS = [80.0, 120.0, 160.0, 200.0]


def nearest(values, target):
    return min(values, key=lambda v: abs(v - target))


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'loopresults.txt'
    df = pd.read_csv(path, sep=r'\s+', skiprows=2, names=COLNAMES)

    lon_values = sorted(df['Lon'].unique())
    f107_values = sorted(df['F107'].unique())
    panels = [nearest(f107_values, f) for f in F107_PANELS]

    # Restrict to the 12 mid-month samples and take the peak over LT
    monthly_doy = (df[['DOY', 'Month']].drop_duplicates()
                     .groupby('Month')['DOY'].first())
    # keep the one DOY per month closest to the 15th, in case the file
    # has more than one sample per month
    peak = (df.groupby(['DOY', 'Lon', 'F107'])['Prob'].max()
              .reset_index())
    peak = peak.merge(df[['DOY', 'Month']].drop_duplicates(), on='DOY')
    monthly = peak.groupby(['Month', 'Lon', 'F107'])['Prob'].mean().reset_index()

    fig, axes = plt.subplots(2, 2, figsize=(18, 14))
    axes = axes.flatten()
    levels = np.linspace(0, 0.8, 17)
    panel_labels = ['(a)', '(b)', '(c)', '(d)']

    cf = None
    for ax, label, f107 in zip(axes, panel_labels, panels):
        sub = monthly[np.isclose(monthly['F107'], f107)]
        pivot = sub.pivot_table(index='Month', columns='Lon',
                                 values='Prob').reindex(
                                     index=range(1, 13), columns=lon_values)

        cf = ax.contourf(lon_values, range(1, 13), pivot.values,
                          levels=levels, cmap='inferno_r', extend='max')
        ax.set_yticks(range(1, 13))
        ax.set_yticklabels(MONTH_NAMES)
        ax.set_xlabel('Longitude (deg.)')
        ax.set_ylabel('Month')
        ax.set_title(f'{label} Monthly IBP index with F10.7 c.m. = '
                      f'{f107:.0f} s.f.u.', fontsize=12)

        # Secondary top axis with sector-name labels
        top_ax = ax.secondary_xaxis('top')
        top_ax.set_xticks(list(SECTOR_LABELS.keys()))
        top_ax.set_xticklabels(list(SECTOR_LABELS.values()), fontsize=9)

        cbar = fig.colorbar(cf, ax=ax, ticks=np.arange(0, 0.9, 0.1))
        cbar.set_label('Occurrence Rate')

    fig.suptitle('IBP Model: Occurrence Probability vs. Month and '
                 'Longitude', fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig('figures/Figure6_F107_Month_Longitude.png',
                dpi=150, bbox_inches='tight')
    print('Wrote figures/Figure6_F107_Month_Longitude.png')


if __name__ == '__main__':
    main()
