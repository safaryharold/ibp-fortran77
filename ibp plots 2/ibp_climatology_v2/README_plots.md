# IBP Climatology Plots — from your `loopresults.txt`

Built directly from the batch run you uploaded: 139,200 rows — a full
factorial of 29 days-of-year (roughly biweekly, with extra points
right at the equinoxes/solstices), 25 longitudes (15° steps, -180°
to +180°), 8 F10.7 levels (80-150 sfu in 10 sfu steps), and 24 local
times (18:00-06:00 next day, 0.5h steps).

## Figures (in `figures/`)

| File | Shows |
|---|---|
| `lt_month_African.png`  | LT vs Month, African sector (0°E) |
| `lt_month_American.png` | LT vs Month, American sector (75°W) |
| `lt_month_Asian.png`    | LT vs Month, Asian sector (90°E) |
| `lt_month_Pacific.png`  | LT vs Month, Pacific sector (165°W) |
| `lt_month_Custom.png`   | LT vs Month, 45°E |
| `lon_month_climatology.png` | Longitude vs Month, all 25 longitudes, fixed LT = 21h |
| `seasonal_profiles.png` | Peak nighttime probability vs month, one panel per sector |
| `f107_response.png`     | **New** — probability vs F10.7 "dose-response" curves, one panel per sector, Jan/Apr/Jul/Oct as separate lines. Only possible because this run has 8 F10.7 levels instead of 3. |

Each heatmap panels three representative F10.7 levels (low/mid/high
from what's actually in your file — currently 80/110/150 sfu) side by
side. The five "sector" longitudes are the closest available values in
your grid to the standard sector centers used in earlier plots.

## Notes on the data

- Multiple days-of-year fall within some months (e.g. 3 samples in
  January, 4 in March); the monthly heatmaps and profiles average
  over whatever days fall in each month. Minor month-to-month texture
  you might notice (e.g. a faint band within June) reflects real
  variation between the sampled days within that month, not a
  plotting artifact.
- Representative longitudes, F10.7 panels, and month labels are all
  auto-detected from your file — this script doesn't assume any
  particular grid, so it'll also work on a different batch run with
  different coverage.

## Regenerating

```bash
python3 plot_from_results.py loopresults.txt
```

Rerun any time with a different results file as the argument.
