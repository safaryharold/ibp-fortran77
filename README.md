# IBP model — FORTRAN 77 implementation

A portable, fixed-form FORTRAN 77 implementation of the **Ionospheric Bubble
Probability (IBP) model** of Stolle et al. (2024, *Space Weather* 22,
e2023SW003809), which gives the probability of post-sunset equatorial plasma
depletions as a function of day of year, longitude, local time and F10.7.

- Compiles with `gfortran -std=legacy`; no external libraries.
- The probability is evaluated analytically through the error function.
- Coefficients (official IBP release v2.0.2) are read at run time from one
  ASCII file, `ibp_coefficients.txt`.
- Agrees with the official IBP v2.0.2 Python package to a maximum absolute
  probability difference of 2.2×10⁻⁶ over 50,000 random inputs (8.2×10⁻⁷ on
  the 504-case test set), entirely due to the error-function approximation.

## Contents

| Folder | Purpose |
|---|---|
| `ibp_single_model_share/` | Library interface (`IBP_INIT`, `IBP_PROBABILITY`) for embedding in a host model, an interactive program and a self-checking test. Start here. |
| `ibp_model_f77/` | Standalone single-run (`pbsingle`) and batch (`pbbatch`) programs. |
| `validation/` | Comparison against the official Python package, an independent Python reference, the coefficient converter and the scripts that produce the paper's figures. See `validation/README.txt`. |

## Quick start

```sh
cd ibp_single_model_share
make            # on macOS with Anaconda: make FC="gfortran -B/usr/bin/"
./test_ibp
```

## Reference

Stolle, C., et al., 2024. Space Weather 22, e2023SW003809.
The accompanying paper describing this implementation is under review.
