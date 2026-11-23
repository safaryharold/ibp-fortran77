"""
plot_from_results.py

Builds climatology/seasonal plots directly from a pbbatch loopresults.txt
file, using whatever grid the file actually contains (auto-detected -
no assumptions about which sectors, F10.7 levels, or LT range were run).

Produces (in figures/):
  1. lt_month_<sector>.png   - LT vs Month heatmap, one per representative
                                longitude, 3 F10.7 panels
  2. lon_month_climatology.png - Longitude vs Month heatmap at a fixed LT,
                                3 F10.7 panels
  3. seasonal_profiles.png   - peak-LT probability vs month, one panel per
                                representative longitude, F10.7 as lines
  4. f107_response.png       - probability vs F10.7 "dose-response" curves,
                                one panel per representative longitude,
                                representative months as lines (only made
                                if the file has more than 3 F10.7 levels -
                                that's what makes this plot worthwhile)

Run with: python3 plot_from_results.py loopresults.txt
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

COLNAMES = ['DOY', 'Lon', 'F107', 'LT', 'Month', 'Prob']
MONTH_NAMES = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
               'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']


def load_results(path):
    df = pd.read_csv(path, sep=r'\s+', skiprows=2, names=COLNAMES)
    return df


def nearest(values, target):
    return min(values, key=lambda v: abs(v - target))


def lt_display_label(lt):
    v = lt % 24
    return f"{int(round(v)):02d}"


def pick_sectors(lon_values):
    """Pick up to 5 representative longitudes closest to the standard
    sector centers actually used in the earlier hand-built grids."""
    targets = {'African': 0.0, 'American': -75.0, 'Asian': 90.0,
               'Pacific': -165.0, 'Custom': 45.0}
    picked = {}
    for name, t in targets.items():
        picked[name] = nearest(lon_values, t)
    return picked


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else 'loopresults.txt'
    df = load_results(path)

    lon_values = sorted(df['Lon'].unique())
    f107_values = sorted(df['F107'].unique())
    lt_values = sorted(df['LT'].unique())
    months_present = sorted(df['Month'].unique())

    print(f"Loaded {len(df)} rows: {len(lon_values)} longitudes, "
          f"{len(f107_values)} F10.7 levels, {len(lt_values)} local "
          f"times, {len(months_present)} months.")

    # Representative F10.7 panels: low / mid / high from what's available
    f107_panels = [f107_values[0],
                    nearest(f107_values, np.median(f107_values)),
                    f107_values[-1]]
    f107_panels = sorted(set(f107_panels))
    while len(f107_panels) < 3 and len(f107_values) >= 3:
        # degenerate case guard; shouldn't trigger for well-formed grids
        break
    f107_labels = {f: f"{f:.0f} sfu" for f in f107_panels}

    sectors = pick_sectors(lon_values)
    print("Representative longitudes:", sectors)
    print("Representative F10.7 panels:", f107_panels)

    # Monthly aggregation: mean over all DOY samples that fall in each
    # calendar month (the grid may sample more than one day per month).
    monthly = (df.groupby(['Month', 'Lon', 'F107', 'LT'])['Prob']
                 .mean().reset_index())

    # ----------------------------------------------------------------
    # 1. LT vs Month heatmaps, one per representative longitude
    # ----------------------------------------------------------------
    for sector_name, lon_val in sectors.items():
        fig, axes = plt.subplots(1, len(f107_panels), figsize=(5*len(f107_panels), 5),
                                  sharey=True)
        if len(f107_panels) == 1:
            axes = [axes]
        sub = monthly[np.isclose(monthly['Lon'], lon_val)]

        im = None
        for ax, f107 in zip(axes, f107_panels):
            s = sub[np.isclose(sub['F107'], f107)]
            pivot = s.pivot_table(index='LT', columns='Month', values='Prob')
            pivot = pivot.reindex(index=lt_values, columns=months_present)

            im = ax.pcolormesh(range(len(months_present) + 1),
                                list(lt_values) + [lt_values[-1] +
                                    (lt_values[1]-lt_values[0])],
                                pivot.values, cmap='viridis',
                                vmin=0, vmax=1, shading='flat')
            ax.set_xticks(np.arange(len(months_present)) + 0.5)
            ax.set_xticklabels([MONTH_NAMES[m-1] for m in months_present],
                                rotation=45, ha='right')
            yticks = [lt for lt in lt_values if lt % 2 == 0]
            ax.set_yticks(yticks)
            ax.set_yticklabels([lt_display_label(lt) for lt in yticks])
            ax.set_title(f107_labels[f107])
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

    # ----------------------------------------------------------------
    # 2. Longitude vs Month heatmap at a fixed representative LT
    # ----------------------------------------------------------------
    lt_fixed = nearest(lt_values, 21.0)
    fig, axes = plt.subplots(1, len(f107_panels), figsize=(5*len(f107_panels), 6),
                              sharey=True)
    if len(f107_panels) == 1:
        axes = [axes]
    sub = monthly[np.isclose(monthly['LT'], lt_fixed)]

    im = None
    for ax, f107 in zip(axes, f107_panels):
        s = sub[np.isclose(sub['F107'], f107)]
        pivot = s.pivot_table(index='Lon', columns='Month', values='Prob')
        pivot = pivot.reindex(index=lon_values, columns=months_present)

        lon_step = lon_values[1] - lon_values[0]
        im = ax.pcolormesh(range(len(months_present) + 1),
                            list(lon_values) + [lon_values[-1] + lon_step],
                            pivot.values, cmap='viridis',
                            vmin=0, vmax=1, shading='flat')
        ax.set_xticks(np.arange(len(months_present)) + 0.5)
        ax.set_xticklabels([MONTH_NAMES[m-1] for m in months_present],
                            rotation=45, ha='right')
        ax.set_yticks(range(-180, 181, 60))
        ax.set_title(f107_labels[f107])
        if ax is axes[0]:
            ax.set_ylabel('Longitude (\u00b0E)')

    fig.colorbar(im, ax=axes, label='Occurrence probability',
                 fraction=0.025, pad=0.02)
    fig.suptitle(f'IBP Occurrence Climatology - Longitude vs Month '
                 f'(fixed LT = {lt_display_label(lt_fixed)}h)', fontsize=14)
    fig.savefig('figures/lon_month_climatology.png',
                dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('Wrote figures/lon_month_climatology.png')

    # ----------------------------------------------------------------
    # 3. Seasonal profile: peak-LT probability vs month, per sector
    # ----------------------------------------------------------------
    n = len(sectors)
    ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.3*ncols, 4.3*nrows))
    axes = np.array(axes).flatten()

    for ax, (sector_name, lon_val) in zip(axes, sectors.items()):
        sub = monthly[np.isclose(monthly['Lon'], lon_val)]
        for f107 in f107_panels:
            s = sub[np.isclose(sub['F107'], f107)]
            peak = s.groupby('Month')['Prob'].max().reindex(months_present)
            ax.plot([MONTH_NAMES[m-1] for m in months_present], peak.values,
                    marker='o', label=f107_labels[f107])
        ax.set_title(f'{sector_name} ({lon_val:+.0f}\u00b0E)')
        ax.set_ylim(0, 1)
        ax.tick_params(axis='x', rotation=45)
        ax.grid(alpha=0.3)

    for ax in axes[:n]:
        pass
    axes[0].set_ylabel('Peak nighttime probability')
    for i in range(n, len(axes)):
        axes[i].axis('off')
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='lower right',
               bbox_to_anchor=(0.98, 0.02))
    fig.suptitle('IBP Seasonal Profile - Peak Nighttime Occurrence '
                 'Probability by Sector', fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig('figures/seasonal_profiles.png', dpi=150, bbox_inches='tight')
    plt.close(fig)
    print('Wrote figures/seasonal_profiles.png')

    # ----------------------------------------------------------------
    # 4. F10.7 dose-response curves (only worthwhile with >3 F10.7 levels)
    # ----------------------------------------------------------------
    if len(f107_values) > 3:
        rep_months = []
        for target in [1, 4, 7, 10]:  # roughly Jan, Apr, Jul, Oct
            rep_months.append(nearest(months_present, target))
        rep_months = sorted(set(rep_months))

        fig, axes = plt.subplots(nrows, ncols, figsize=(5.3*ncols, 4.3*nrows))
        axes = np.array(axes).flatten()

        for ax, (sector_name, lon_val) in zip(axes, sectors.items()):
            sub = monthly[np.isclose(monthly['Lon'], lon_val)]
            for m in rep_months:
                s = sub[sub['Month'] == m]
                peak = s.groupby('F107')['Prob'].max().reindex(f107_values)
                ax.plot(f107_values, peak.values, marker='o',
                        label=MONTH_NAMES[m-1])
            ax.set_title(f'{sector_name} ({lon_val:+.0f}\u00b0E)')
            ax.set_ylim(0, 1)
            ax.set_xlabel('F10.7 (sfu)')
            ax.grid(alpha=0.3)

        axes[0].set_ylabel('Peak nighttime probability')
        for i in range(n, len(axes)):
            axes[i].axis('off')
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(handles, labels, loc='lower right',
                   bbox_to_anchor=(0.98, 0.02), title='Month')
        fig.suptitle('IBP Solar-Flux Response - Peak Nighttime Probability '
                     'vs F10.7 by Sector', fontsize=14)
        fig.tight_layout(rect=[0, 0, 1, 0.96])
        fig.savefig('figures/f107_response.png', dpi=150, bbox_inches='tight')
        plt.close(fig)
        print('Wrote figures/f107_response.png')
    else:
        print('Skipping F10.7 dose-response plot (only',
              len(f107_values), 'F10.7 levels in this file - not worth it).')

    print('\nAll figures written to figures/')


if __name__ == '__main__':
    main()
