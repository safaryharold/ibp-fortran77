# IBP Single-Run Model — Sharing Package

The Ionospheric Bubble Probability (IBP) single-run model, packaged
for handing off to another team / incorporating into a larger model.
FORTRAN 77 (fixed-form), no external dependencies beyond a Fortran
compiler.

## What's in here

| File | What it is |
|---|---|
| `ibp_probability.f` | **The library interface** — `IBP_INIT` + `IBP_PROBABILITY`. No stdout/stderr output or `STOP`; initialization reports coefficient-file errors through `IER`. |
| `test_ibp_probability.f` | Minimal working example of calling the library, with a self-check against known-good reference values. Read this first. |
| `pbsingle.f` | The interactive, standalone version — prompts for input, prints full debug detail to the screen. Useful for manual testing/demos, not for embedding. |
| `pbutils.f` | Shared physics functions and the coefficient-file loader (`READCF`). Required by both `ibp_probability.f` and `pbsingle.f`. |
| `pbcoeff.inc` | `COMMON` block holding the model coefficients in memory. `INCLUDE`d by the files above. |
| `ibp_coefficients.txt` | The model coefficients themselves, as plain text (not compiled in — loaded at run time). Must ship alongside the compiled model; the path is passed to `IBP_INIT` at run time, so it can live wherever your host model expects. |
| `Makefile` | Builds everything (see below). |

## Quick start (for a Fortran host model)

```fortran
      CALL IBP_INIT('ibp_coefficients.txt', IER)
      IF (IER .NE. 0) STOP 'could not load IBP coefficients'

      ! ... inside whatever loop your model already has ...
      CALL IBP_PROBABILITY(IDOY, ALON, F107IN, TIN,
     +                      PROB, IER, WARN_DOY, WARN_F107, WARN_Y0)
```

`IBP_INIT` is called once, at startup. Check `IER` before calling
`IBP_PROBABILITY`; a failed or omitted initialization leaves the model
coefficients unusable. `IBP_PROBABILITY` does not write output or
change the coefficients, so it can be called repeatedly after a
successful initialization. The coefficient table is stored in a
Fortran `COMMON` block: initialize once before parallel work, and do
not call `IBP_INIT` concurrently with probability calculations.

Build the library and run the example:
```bash
make            # builds pbsingle, libibp.a, and test_ibp
make test       # runs the self-check
```

To link `libibp.a` into your own program: compile your code against
`pbcoeff.inc` if it needs the coefficients directly (it normally
won't — the interface below is all you need), then link with
`-L<path> -libp`.

## Interface

```fortran
      SUBROUTINE IBP_INIT(CFFILE, IER)
```
- `CFFILE` (in, `CHARACTER*(*)`) — path to `ibp_coefficients.txt`.
- `IER` (out, `INTEGER`) — 0 on success; 1 if the file cannot be
  opened or any required record cannot be read. No `STOP` is issued.
  Call once before any call to `IBP_PROBABILITY`.

```fortran
      SUBROUTINE IBP_PROBABILITY(IDOY, ALON, F107IN, TIN,
     +                           PROB, IER, WARN_DOY, WARN_F107, WARN_Y0)
```

| Argument | Type | In/Out | Meaning |
|---|---|---|---|
| `IDOY` | `INTEGER` | in | Day of year, nominally 1-365 |
| `ALON` | `DOUBLE PRECISION` | in | Longitude, degrees East (any value; wraps automatically) |
| `F107IN` | `DOUBLE PRECISION` | in | F10.7 solar flux, sfu (nominal valid range 80-200) |
| `TIN` | `DOUBLE PRECISION` | in | Local time, hours (any value; wraps automatically, e.g. 25.5 = 01:30) |
| `PROB` | `DOUBLE PRECISION` | out | Occurrence probability, in [0,1] |
| `IER` | `INTEGER` | out | 0 = success. 1 = coefficients not loaded/initialized or invalid `TAU` (mean lifetime, = 1/γ) — `PROB` is 0.0 and **not physically meaningful**. 2 = non-finite input — `PROB` is 0.0. |
| `WARN_DOY` | `LOGICAL` | out | `.TRUE.` if `IDOY` was outside 1-365; defaulted to January and proceeded |
| `WARN_F107` | `LOGICAL` | out | `.TRUE.` if `F107IN` was outside 80-200 sfu; proceeded anyway |
| `WARN_Y0` | `LOGICAL` | out | `.TRUE.` if an internal term was clamped to avoid overflow (see below) |

Neither public routine writes to stdout/stderr or calls `STOP`.
Initialization file errors and calculation status are returned through
`IER`; input range conditions use the `WARN_*` flags.

## Two deliberate differences from `pbsingle.f`

`ibp_probability.f` produces the same numbers as `pbsingle.f` for
in-range inputs, but isn't a line-for-line port. Two things were
changed specifically because this is meant to run *inside* another
model rather than as a standalone diagnostic tool:

1. **No `STOP`.** `pbsingle.f` calls `STOP` if a coefficient (`TAU`,
   the mean bubble lifetime 1/γ) is invalid. That's fine for a one-off interactive run; it would be
   unacceptable in a library routine, since it would halt the entire
   host model over one bad calculation. This is now `IER=1` instead.

2. **Overflow is clamped, and flagged rather than printed.** Deep in
   the probability calculation there's a term `EXP(y0)` that can
   overflow for extreme inputs. All three drivers (`pbbatch.f`,
   `pbsingle.f`, and this library routine) clamp it the same way; the
   only difference here is that the library raises `WARN_Y0` instead
   of printing to stdout, so a host model can log that it happened.

3. **Initialization is verified.** `IBP_PROBABILITY` checks the
   `IBP_INITED` flag set by a successful `IBP_INIT`; calling it
   before a successful initialization returns `IER=1` instead of
   using whatever happens to be in the `COMMON` block.

4. **Non-finite inputs are rejected.** NaN or non-finite
   `ALON`/`F107IN`/`TIN` return `IER=2` with `PROB=0.0` rather than
   propagating through the trig/exponential terms.

`READCF` also picks the first unopened unit in 90-99 (instead of a
hard-coded unit 12 that could collide with host I/O) and rejects
physically impossible coefficients — non-positive `TAU`/`SIGMA` or
negative/non-finite density-table entries — before marking the
common block usable.

A smaller, related fix: when `IDOY` is out of range, `pbsingle.f`
defaults the *month* to January but still feeds the original,
invalid `IDOY` into the seasonal-oscillation calculation — which both
extrapolates a trig function outside its intended domain and prints
its own separate warning (a second hidden I/O side effect). Here,
an out-of-range `IDOY` is sanitized to 1 everywhere, consistently,
before it's used for anything, and no output is ever printed.

## Validation

`test_ibp_probability.f` checks that a missing coefficient file is
reported without stopping, then checks 9 reference cases (spanning different
seasons, hemispheres, longitude wraparound, solar activity levels,
and month-boundary DOYs 60/91/335 plus out-of-range DOY 366) against
expected probabilities — matches are within 1×10⁻⁵ in every case, and
the out-of-range cases confirm `WARN_DOY` fires without halting the
program. Run `make test` to reproduce this yourself.

## Notes

Coefficients, including the 360×12 longitude/month table, are read
from `ibp_coefficients.txt` (generated from the official IBP v2.0.2
coefficient CDF by `validation/cdf_to_coefficients.py`). Residual
differences from the official Python package (on the order of 10⁻⁶ in
probability) come from `ERF77`, an Abramowitz & Stegun approximation
used because plain FORTRAN 77 has no `ERF` intrinsic (it was only
added in Fortran 2008). Day of year is mapped to month with the same
non-leap calendar as the official package.
