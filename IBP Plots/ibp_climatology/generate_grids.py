"""
generate_grids.py

Builds the two batch-input grids used for the IBP climatology plots.
Edit the constants below to change coverage, then re-run this script
followed by pbbatch (once per grid) and plot_climatology.py.

Outputs:
  batch_inputs_sector.txt  - for the LT-vs-Month / per-sector plots
  batch_inputs_lonmap.txt  - for the Longitude-vs-Month plot
  meta.json                 - shared metadata read by plot_climatology.py
"""

import datetime
import json

# ------------------------------------------------------------------
# Edit these to change coverage
# ------------------------------------------------------------------
YEAR = 2025  # non-leap reference year, only used to compute mid-month DOY

F107_LEVELS = [70.0, 120.0, 180.0]  # low / moderate / high solar flux

SECTORS = {
    'African':  0.0,
    'American': -75.0,
    'Asian':    100.0,
    'Pacific':  -160.0,
    'Custom45': 45.0,
}

LT_START, LT_END, LT_STEP = 18.0, 30.0, 0.5   # local time range for sector grid
LON_START, LON_END, LON_STEP = -180, 180, 5   # longitude range for lonmap grid
LT_FIXED = 21.0                                # fixed LT for the lonmap grid

# ------------------------------------------------------------------
# Build mid-month day-of-year list
# ------------------------------------------------------------------
mid_doy, month_names = [], []
for m in range(1, 13):
    d = datetime.date(YEAR, m, 15)
    mid_doy.append(d.timetuple().tm_yday)
    month_names.append(d.strftime('%b'))

# ------------------------------------------------------------------
# Grid A: LT vs Month, per longitude sector
# ------------------------------------------------------------------
lt_values = []
lt = LT_START
while lt <= LT_END + 1e-9:
    lt_values.append(round(lt, 2))
    lt += LT_STEP

rows_a = []
for doy in mid_doy:
    for lon in SECTORS.values():
        for f107 in F107_LEVELS:
            for ltv in lt_values:
                rows_a.append(f"{doy} {lon} {f107} {ltv}")

with open('batch_inputs_sector.txt', 'w') as f:
    f.write('\n'.join(rows_a) + '\n')
print(f"Wrote batch_inputs_sector.txt ({len(rows_a)} rows)")

# ------------------------------------------------------------------
# Grid B: Longitude vs Month, fixed LT
# ------------------------------------------------------------------
lon_values = list(range(LON_START, LON_END + 1, LON_STEP))

rows_b = []
for doy in mid_doy:
    for lon in lon_values:
        for f107 in F107_LEVELS:
            rows_b.append(f"{doy} {float(lon)} {f107} {LT_FIXED}")

with open('batch_inputs_lonmap.txt', 'w') as f:
    f.write('\n'.join(rows_b) + '\n')
print(f"Wrote batch_inputs_lonmap.txt ({len(rows_b)} rows)")

# ------------------------------------------------------------------
# Metadata for plot_climatology.py
# ------------------------------------------------------------------
meta = {
    'month_names': month_names,
    'mid_doy': mid_doy,
    'f107_levels': F107_LEVELS,
    'sectors': SECTORS,
    'lt_values': lt_values,
    'lon_values': [float(x) for x in lon_values],
    'lt_fixed': LT_FIXED,
}
with open('meta.json', 'w') as f:
    json.dump(meta, f, indent=2)
print("Wrote meta.json")
