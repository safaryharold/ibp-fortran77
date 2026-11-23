"""
plot_lt_vs_longitude_seasons.py

Reproduces the "IBP Occurrence Probability: Local Time vs. Longitude"
style figure: a 2x2 grid of filled-contour panels (December Solstice,
March Equinox, June Solstice, September Equinox) at a fixed F10.7,
with Longitude on the x-axis and Local Time on the y-axis.

Run with: python3 plot_lt_vs_longitude_seasons.py [loopresults.txt] [F10.7]
Requires: pandas, numpy, matplotlib
"""

import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

COLNAMES = ['DOY', 'Lon', 'F107', 'LT', 'Month', 'Prob']

# Which season maps to which DOY - edit if you used different reference
# dates when generating your own grid.
SEASONS = {
    'December Solstice': 355,
    'March Equinox': 79,
    'June Solstice': 172,
    'September Equinox': 265,
}
PANEL_ORDER = ['December Solstice', 'March Equinox',
               'June Solstice', 'September Equinox']
PANEL_LABELS = ['(a)', '(b)', '(c)', '(d)']


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'loopresults.txt'
    f107_target = float(sys.argv[2]) if len(sys.argv) > 2 else 150.0

    df = pd.read_csv(path, sep=r'\s+', skiprows=2, names=COLNAMES)

    lon_values = sorted(df['Lon'].unique())
    lt_values = sorted(df['LT'].unique())
    f107_actual = min(df['F107'].unique(), key=lambda v: abs(v - f107_target))

    fig, axes = plt.subplots(2, 2, figsize=(16, 14))
    axes = axes.flatten()
    levels = np.linspace(0, 0.9, 21)

    cf = None
    for ax, label, season_name in zip(axes, PANEL_LABELS, PANEL_ORDER):
        doy = SEASONS[season_name]
        sub = df[np.isclose(df['F107'], f107_actual) & (df['DOY'] == doy)]
        pivot = sub.pivot_table(index='LT', columns='Lon',
                                 values='Prob').reindex(index=lt_values,
                                                          columns=lon_values)

        cf = ax.contourf(lon_values, lt_values, pivot.values,
                          levels=levels, cmap='RdBu_r', extend='max')
        ax.contour(lon_values, lt_values, pivot.values,
                   levels=levels, colors='white', linewidths=0.3, alpha=0.5)
        ax.set_title(f'{label} {season_name}', fontsize=13)
        ax.set_xlabel('Longitude (deg.)')
        ax.set_ylabel('Local Time (h)')

        cbar = fig.colorbar(cf, ax=ax, ticks=np.arange(0, 1.0, 0.15))
        cbar.set_label('Occurrence Rate')

    fig.suptitle('IBP Occurrence Probability: Local Time vs. Longitude\n'
                 f'F10.7 = {f107_actual:.0f} s.f.u', fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    fig.savefig('figures/Figure5_Season_LT_Longitude.png',
                dpi=150, bbox_inches='tight')
    print('Wrote figures/Figure5_Season_LT_Longitude.png')


if __name__ == '__main__':
    main()
