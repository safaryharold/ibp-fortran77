"""
plot_climatology.py

Reads the two batch-mode result files produced by pbbatch:
  - loopresults_sector.txt  (LT x Month, per longitude sector, 3 F10.7 levels)
  - loopresults_lonmap.txt  (Longitude x Month at a fixed LT, 3 F10.7 levels)

and produces:
  1. One LT-vs-Month occurrence-probability heatmap per longitude/sector
     (5 figures: African, American, Asian, Pacific, Custom45), each with
     three side-by-side panels for low/moderate/high F10.7.
  2. One Longitude-vs-Month heatmap at the fixed reference local time,
     with the same three F10.7 panels.
  3. One seasonal-profile line plot (peak nighttime probability vs month)
     with one panel per longitude/sector, each showing all three F10.7
     levels as separate lines.

Run with: python3 plot_climatology.py
Requires: pandas, numpy, matplotlib
"""

import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ------------------------------------------------------------------
# Load metadata and results
# ------------------------------------------------------------------
with open('meta.json') as f:
    META = json.load(f)

MONTH_NAMES = META['month_names']
MID_DOY     = META['mid_doy']
F107_LEVELS = META['f107_levels']
SECTORS     = META['sectors']        # name -> longitude
LT_VALUES   = META['lt_values']
LON_VALUES  = META['lon_values']
LT_FIXED    = META['lt_fixed']

DOY_TO_MONTH = dict(zip(MID_DOY, MONTH_NAMES))
F107_LABELS = {70.0: 'Low (70 sfu)', 120.0: 'Moderate (120 sfu)',
               180.0: 'High (180 sfu)'}

COLNAMES = ['DOY', 'Lon', 'F107', 'LT', 'Month', 'Prob']


def load_results(path):
    df = pd.read_csv(path, sep=r'\s+', skiprows=2, names=COLNAMES)
    df['MonthName'] = df['DOY'].map(DOY_TO_MONTH)
    return df


sector_df = load_results('loopresults_sector.txt')
lonmap_df = load_results('loopresults_lonmap.txt')


def lt_display_label(lt):
    """Map an internal LT value (18..30) to a 0-23h clock label."""
    v = lt % 24
    return f"{int(round(v)):02d}"


# ------------------------------------------------------------------
# 1. LT vs Month heatmaps, one figure per sector/longitude
# ------------------------------------------------------------------
for sector_name, lon_val in SECTORS.items():
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    sub = sector_df[np.isclose(sector_df['Lon'], lon_val)]

    for ax, f107 in zip(axes, F107_LEVELS):
        s = sub[np.isclose(sub['F107'], f107)]
        pivot = s.pivot_table(index='LT', columns='DOY', values='Prob')
        pivot = pivot.reindex(index=LT_VALUES, columns=MID_DOY)

        im = ax.pcolormesh(range(len(MID_DOY) + 1),
                            list(LT_VALUES) + [LT_VALUES[-1] + 0.5],
                            pivot.values, cmap='viridis',
                            vmin=0, vmax=1, shading='flat')
        ax.set_xticks(np.arange(len(MID_DOY)) + 0.5)
        ax.set_xticklabels(MONTH_NAMES, rotation=45, ha='right')
        yticks = [lt for lt in LT_VALUES if lt % 2 == 0]
        ax.set_yticks(yticks)
        ax.set_yticklabels([lt_display_label(lt) for lt in yticks])
        ax.set_title(F107_LABELS[f107])
        if ax is axes[0]:
            ax.set_ylabel('Local Time (h)')

    fig.colorbar(im, ax=axes, label='Occurrence probability',
                 fraction=0.025, pad=0.02)
    fig.suptitle(f'IBP Occurrence Climatology - {sector_name} sector '
                 f'({lon_val:+.0f}\u00b0E)', fontsize=14)
    fig.savefig(f'figures/lt_month_{sector_name}.png',
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f'Wrote figures/lt_month_{sector_name}.png')


# ------------------------------------------------------------------
# 2. Longitude vs Month heatmap at fixed LT
# ------------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)

for ax, f107 in zip(axes, F107_LEVELS):
    s = lonmap_df[np.isclose(lonmap_df['F107'], f107)]
    pivot = s.pivot_table(index='Lon', columns='DOY', values='Prob')
    pivot = pivot.reindex(index=LON_VALUES, columns=MID_DOY)

    im = ax.pcolormesh(range(len(MID_DOY) + 1),
                        list(LON_VALUES) + [LON_VALUES[-1] + 5],
                        pivot.values, cmap='viridis',
                        vmin=0, vmax=1, shading='flat')
    ax.set_xticks(np.arange(len(MID_DOY)) + 0.5)
    ax.set_xticklabels(MONTH_NAMES, rotation=45, ha='right')
    ax.set_yticks(range(-180, 181, 60))
    ax.set_title(F107_LABELS[f107])
    if ax is axes[0]:
        ax.set_ylabel('Longitude (\u00b0E)')

fig.colorbar(im, ax=axes, label='Occurrence probability',
             fraction=0.025, pad=0.02)
fig.suptitle(f'IBP Occurrence Climatology - Longitude vs Month '
             f'(fixed LT = {LT_FIXED:.0f}h)', fontsize=14)
fig.savefig('figures/lon_month_climatology.png',
            dpi=150, bbox_inches='tight')
plt.close(fig)
print('Wrote figures/lon_month_climatology.png')


# ------------------------------------------------------------------
# 3. Seasonal profile: peak-LT probability vs month, per sector
# ------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(16, 9))
axes = axes.flatten()

for ax, (sector_name, lon_val) in zip(axes, SECTORS.items()):
    sub = sector_df[np.isclose(sector_df['Lon'], lon_val)]
    for f107 in F107_LEVELS:
        s = sub[np.isclose(sub['F107'], f107)]
        peak = s.groupby('DOY')['Prob'].max().reindex(MID_DOY)
        ax.plot(MONTH_NAMES, peak.values, marker='o',
                label=F107_LABELS[f107])
    ax.set_title(f'{sector_name} ({lon_val:+.0f}\u00b0E)')
    ax.set_ylim(0, 1)
    ax.tick_params(axis='x', rotation=45)
    ax.grid(alpha=0.3)

axes[0].set_ylabel('Peak nighttime probability')
axes[3].set_ylabel('Peak nighttime probability')
axes[-1].axis('off')
handles, labels = axes[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='lower right',
           bbox_to_anchor=(0.98, 0.05))
fig.suptitle('IBP Seasonal Profile - Peak Nighttime Occurrence '
             'Probability by Sector', fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig('figures/seasonal_profiles.png', dpi=150, bbox_inches='tight')
plt.close(fig)
print('Wrote figures/seasonal_profiles.png')

print('\nAll figures written to figures/')
