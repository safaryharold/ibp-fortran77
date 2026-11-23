"""
generate_grid.py

Builds the single comprehensive batch-input grid used by all three
plotting scripts in this package. Edit the constants below to change
coverage (finer longitude step, different F10.7 levels, etc.), then
re-run this script followed by ./pbbatch and each plot_*.py script.

Output: batch_inputs.txt, meta.json
"""

import datetime
import json

# ------------------------------------------------------------------
# Edit these to change coverage
# ------------------------------------------------------------------
YEAR = 2025  # non-leap reference year, only used to compute DOYs

LON_START, LON_END, LON_STEP = -180, 180, 5     # degrees
LT_START, LT_END, LT_STEP = 18.0, 29.5, 0.5      # hours (18:00-06:00 next day)
F107_VALUES = [80.0, 120.0, 150.0, 160.0, 200.0]  # sfu

# Reference dates for the seasonal (Figure 5 style) panels
SEASON_DATES = {
    'December Solstice': (12, 21),
    'March Equinox':     (3, 20),
    'June Solstice':     (6, 21),
    'September Equinox': (9, 22),
}

# ------------------------------------------------------------------
# Build DOY lists
# ------------------------------------------------------------------
mid_doy, month_names = [], []
for m in range(1, 13):
    d = datetime.date(YEAR, m, 15)
    mid_doy.append(d.timetuple().tm_yday)
    month_names.append(d.strftime('%b'))

seasons = {}
for name, (m, d) in SEASON_DATES.items():
    seasons[name] = datetime.date(YEAR, m, d).timetuple().tm_yday

all_doy = sorted(set(mid_doy) | set(seasons.values()))

lon_values = [float(x) for x in range(LON_START, LON_END + 1, LON_STEP)]
lt_values = []
lt = LT_START
while lt <= LT_END + 1e-9:
    lt_values.append(round(lt, 2))
    lt += LT_STEP

# ------------------------------------------------------------------
# Write batch_inputs.txt
# ------------------------------------------------------------------
rows = []
for doy in all_doy:
    for lon in lon_values:
        for f107 in F107_VALUES:
            for ltv in lt_values:
                rows.append(f"{doy} {lon} {f107} {ltv}")

with open('batch_inputs.txt', 'w') as f:
    f.write('\n'.join(rows) + '\n')
print(f"Wrote batch_inputs.txt ({len(rows)} rows)")

meta = {
    'month_names': month_names, 'mid_doy': mid_doy,
    'seasons': seasons, 'all_doy': all_doy,
    'lon_values': lon_values, 'lt_values': lt_values,
    'f107_values': F107_VALUES,
}
with open('meta.json', 'w') as f:
    json.dump(meta, f, indent=2)
print("Wrote meta.json")
