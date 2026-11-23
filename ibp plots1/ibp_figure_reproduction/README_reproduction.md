# Reproducing Figures 3, 5, and 6 from the F77 IBP Model

This package rebuilds the style of three reference figures you
shared, driven by the F77 `pbbatch` model instead of whatever
originally produced them. Since I don't have your original plotting
code, a few things had to be *decided* rather than copied exactly —
they're called out below so you know what's a faithful reproduction
and what's my interpretation.

## Files

| File | Purpose |
|---|---|
| `generate_grid.py` | Builds `batch_inputs.txt` (140,160 rows) covering everything all three figures need |
| `batch_inputs.txt` / `loopresults.txt` | The grid and the model's output on it |
| `plot_lt_vs_longitude_seasons.py` | → `Figure5_Season_LT_Longitude.png` |
| `plot_month_vs_longitude_f107panels.py` | → `Figure6_F107_Month_Longitude.png` |
| `plot_monthly_longitude_density.py` | → `Figure3_Monthly_Longitude_Density.png` |
| `pbbatch`, `ibp_coefficients.txt` | Copied in so the package is self-contained |

## The grid

One grid feeds all three figures:
- **Longitude**: -180° to 180°, 5° steps (73 points) — fine enough for
  smooth-looking contours, since the model's own longitude table is
  1°-resolution internally.
- **Local time**: 18:00–06:00 (next day), 0.5h steps (24 points).
- **F10.7**: 80, 120, 150, 160, 200 sfu — covers both Figure 5 (150)
  and Figure 6 (80/120/160/200).
- **Day-of-year**: the 12 mid-month days (for Figures 3 & 6) plus the
  4 exact 2025 solstice/equinox dates (for Figure 5): Mar 20, Jun 21,
  Sep 22, Dec 21.

## Decisions I made (things your original figures don't tell me)

- **Figure 5** (`plot_lt_vs_longitude_seasons.py`): straightforward —
  the axes and F10.7 are exactly what the reference figure shows, so
  this one is a direct reproduction, just at whatever exact
  solstice/equinox date you use.

- **Figure 6** (`plot_month_vs_longitude_f107panels.py`): the
  reference figure's "Monthly IBP index" isn't something the F77
  model outputs directly — the model gives probability *per local
  time*, not a single monthly number. I used **peak probability over
  the local-time range** for each (month, longitude, F10.7) as that
  index. If your original used something else (e.g. an average over
  local time, or a specific LT), edit the `peak = df.groupby(...)`
  line to `.mean()` instead of `.max()`, or filter to a fixed LT
  first. The sector labels (Pacific/America/Atlantic Africa/India)
  are placed at round longitude values (-180°, -90°, 0°, 80°, 180°)
  matching the reference's approximate positions — adjust
  `SECTOR_LABELS` if you know the exact boundaries your original used.

- **Figure 3** (`plot_monthly_longitude_density.py`): same "no
  built-in index" issue, plus this figure needs a **weighting**
  scheme to turn probability into a density. I used peak-over-LT
  probability at a single reference F10.7 (120 sfu, passed as a
  command-line argument — change it freely) as the per-longitude
  *weight*, then normalized so each month's histogram integrates to
  1 and overlaid a weighted Gaussian KDE. If your original weighted
  by something else (e.g. summed over all F10.7 levels, or used raw
  counts rather than probability-as-weight), that's the line to
  change (`monthly = peak.groupby(...)`).

None of these change the underlying model output — only how it's
aggregated for these three specific plots. If you can tell me more
about how the original figures were built (e.g. a script, or even
just "we averaged over LT" / "we used F10.7=100 for figure 3"), I can
tighten these to match exactly instead of approximating.

## Regenerating or customizing

```bash
python3 generate_grid.py        # rewrite batch_inputs.txt / meta.json
./pbbatch                       # rerun the model (rename loopresults.txt if needed)
python3 plot_lt_vs_longitude_seasons.py loopresults.txt 150
python3 plot_month_vs_longitude_f107panels.py loopresults.txt
python3 plot_monthly_longitude_density.py loopresults.txt 120
```

Each plotting script takes the results file (and, where relevant, a
reference F10.7) as command-line arguments, so you can point them at
any other `loopresults.txt` — including the ones from earlier in this
conversation — without editing the scripts.
