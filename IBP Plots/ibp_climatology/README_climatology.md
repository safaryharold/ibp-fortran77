# IBP Climatology & Seasonal Plots

Generated from the F77 `pbbatch` model with three purpose-built input
grids (the 8-row sample file isn't enough for climatology — these runs
use 4,500 + 2,628 rows respectively).

## Figures (in `figures/`)

| File | Shows |
|---|---|
| `lt_month_African.png`  | LT vs Month occurrence probability, African sector (0°E) |
| `lt_month_American.png` | LT vs Month, American sector (75°W) |
| `lt_month_Asian.png`    | LT vs Month, Asian sector (100°E) |
| `lt_month_Pacific.png`  | LT vs Month, Pacific sector (160°W) |
| `lt_month_Custom45.png` | LT vs Month, single longitude example (45°E) |
| `lon_month_climatology.png` | Longitude vs Month, all longitudes, fixed LT = 21h |
| `seasonal_profiles.png` | Peak nighttime probability vs month, one panel per sector |

Every figure shows all three F10.7 levels (70 / 120 / 180 sfu) side by
side (or as separate lines) so solar-activity dependence is directly
comparable.

## How the grids were built

- **Months**: mid-month day-of-year for a non-leap reference year
  (Jan 15 = DOY 15, ... Dec 15 = DOY 349).
- **Sector grid** (`batch_inputs_sector.txt`, → `loopresults_sector.txt`):
  12 months × 25 local times (18.0–30.0h in 0.5h steps, i.e. 18:00
  through 06:00 next day) × 5 longitudes × 3 F10.7 levels = 4,500 runs.
- **Longitude-map grid** (`batch_inputs_lonmap.txt`, →
  `loopresults_lonmap.txt`): 12 months × 73 longitudes (-180° to
  +180° in 5° steps) × 3 F10.7 levels, at a fixed LT = 21:00
  = 2,628 runs.

## Regenerating or customizing

1. `python3 generate_grids.py` — rewrites the two `batch_inputs_*.txt`
   files. Edit the constants at the top (sectors, F10.7 levels, LT
   range, longitude step) to change coverage.
2. Copy each grid to `batch_inputs.txt`, run `./pbbatch` (the compiled
   F77 executable — build it from the `ibp_model_f77` package if you
   don't already have it), then rename `loopresults.txt` to match
   (`loopresults_sector.txt` or `loopresults_lonmap.txt`).
3. `python3 plot_climatology.py` — regenerates all 7 PNGs into
   `figures/`.

To add a sector or move the custom longitude, edit the `SECTORS` dict
in `generate_grids.py`; to change the fixed LT for the longitude map,
edit `LT_FIXED`. Both scripts read/write `meta.json`, which the
plotting script uses to know what was actually run — keep them in
sync (i.e. re-run `plot_climatology.py` any time you regenerate the
grids).
