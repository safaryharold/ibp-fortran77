# IBP Model — FORTRAN 77 Version

This is a fixed-form FORTRAN 77 implementation of the Ionospheric
Bubble Probability (IBP) model (Stolle et al., 2024) for batch and
interactive use. It has been validated against the official IBP
Python package (see "Accuracy and validation" below).

## Files in this package

| File                     | Purpose                                                          |
|--------------------------|-------------------------------------------------------------------|
| `pbbatch.f`              | Batch driver |
| `pbsingle.f`             | Interactive single-run driver |
| `pbutils.f`              | Coefficient loader + physics functions, shared by both programs |
| `pbcoeff.inc`            | F77 `INCLUDE` file — the `COMMON` block holding the coefficients |
| `ibp_coefficients.txt`   | **The model coefficients (IBP v2.0.2), in a plain-text data file read at run time** |
| `batch_inputs.txt`       | Sample input for `pbbatch` (8 rows) — replace with your own data |
| `Makefile`                | Builds everything with one command |
| `README.md`              | This file |

`pbbatch` and `pbsingle` are two separate main programs that both
link against the same `pbutils.f` and share the same
`ibp_coefficients.txt`, so a coefficient update affects both
consistently.

## Two programs, one shared physics/coefficients layer

- **`pbbatch`** — reads rows from `batch_inputs.txt`, writes
  `loopresults.txt`. Auto-detects and skips a header row if present.
- **`pbsingle`** — interactive: prompts you for DOY, longitude,
  F10.7, and local time one at a time, then prints the same
  calculation with full debug/intermediate-value detail to the
  screen (no files read or written besides the coefficients file).
  Like `pbbatch`, this program clamps the `y0` overflow term before
  taking `EXP(y0)` (for `|y0| > 100`) and prints a warning when it
  does so. All three drivers (`pbbatch`, `pbsingle`, and the
  `ibp_probability` library routine) now share the identical clamp,
  so they agree even for extreme synthetic inputs. Note that for all
  physically meaningful inputs `|y0|` stays well below 100, so the
  clamp never actually engages in normal use.

Both share `pbutils.f` (the physics functions and `READCF`
coefficient loader) and `pbcoeff.inc` (the `COMMON` block), and both
read coefficients from the same `ibp_coefficients.txt`.

## FORTRAN 77 implementation notes

FORTRAN 77 has no modules, derived types, allocatable arrays,
internal procedures, `CYCLE`/`EXIT`, `MODULO`, `TRIM` or `ERF`. The
code handles this as follows:

- Coefficients live in a `COMMON /PBCOEF/` block declared once in
  `pbcoeff.inc` and `INCLUDE`d where needed, with a fixed-size
  `PHITAB(360,12)` density table.
- Day of year → month is the external function `MONOF`
  (non-leap calendar).
- Loops use labelled `DO`/`CONTINUE` with `GOTO`.
- `DMODLO` provides `MODULO`, built from the standard `DMOD`.
- `ERF77` provides the error function using the Abramowitz & Stegun
  7.1.26 rational approximation (max absolute error 1.5×10⁻⁷).
- End of file is detected from a negative `IOSTAT`.
- `FSTNB` / `LSTNB` find the first / last non-blank character, used
  to detect and report an optional header row in the input file.

## Header-line auto-detection

`pbbatch.f` peeks at the first line of
`batch_inputs.txt`: if it starts with a digit or a minus sign, it's
treated as data and the file is rewound; otherwise it's treated as a
header row, reported once (`Header line detected and skipped: ...`),
and left consumed. You can hand it a file with or without a header
row and it will do the right thing either way.

## Coefficients: a plain-text file, not source code

All coefficients, including the 360×12 `phi_table`, live in
`ibp_coefficients.txt`, loaded at run time by the `READCF` subroutine
in `pbutils.f`. The file is generated from the official IBP
coefficient CDF (`..._0101.cdf`, IBP v2.0.2) by
`validation/cdf_to_coefficients.py`, at full double precision, with
each value right-justified in a fixed-width field so columns line up.
To use another IBP release, convert its CDF with that script; no
recompilation is needed. The layout is:

```
<title line>
SCALAR PARAMS: C1 C2 TAU(=1/gamma,mean_lifetime_h) MU0 SIGMA
<line of 5 numbers>
DIURNAL OSC PARAMS: G0 GS(1) GC(1) GS(2) GC(2)
<line of 5 numbers>
MONTHLY LT-SHIFT T0MON (Jan..Dec, 12 values)
<line of 12 numbers>
... (4 more 12-value monthly harmonic lines) ...
LON-MONTH PHI TABLE: 360 rows x 12 cols
<360 lines of 12 numbers each>
```

The section-header lines (the plain-text descriptions) are read and
discarded by `READCF` — they're there purely so a human opening the
file can tell what each block is. If you need to update a
coefficient, regenerate or edit `ibp_coefficients.txt`; no
recompilation is needed. If you add or remove a header/comment line, you must add or
remove the matching `READ(IUCF,'(A)') LINE` in `READCF` (in
`pbutils.f`) to keep the file and the reader in sync — F77 has no way
to auto-detect comment lines the way modern parsers can.

## Accuracy and validation

**Primary validation — against the reference model.** The FORTRAN 77
output was validated against the official IBP v2.0.2 Python
implementation over the 504-case grid in `validation/`. The maximum
absolute probability difference is 8.2×10⁻⁷ and the mean is 4.4×10⁻⁸
(this is the number quoted in the paper). It is fully reproducible from
the supplied files: run `python validation/compare_504.py`, which
evaluates an independent Python reference (`validation/ibp_reference.py`,
reading the same `ibp_coefficients.txt`) over the grid and differences
it against `fortran_loopresults_504.txt`. The residual is dominated by
the `ERF77` error-function approximation — the reference uses exact
`erf` — exactly as described above. See `validation/README.txt`.

**Second check — against the official package directly.**
`validation/compare_official.py` (with the `validation/ibp_dump.f`
driver) evaluates 50,000 random inputs with the Fortran library and
with the official `ibpmodel` 2.0.2 code. The maximum difference is
2.2×10⁻⁶ in probability; with `ERF77` substituted into the official
code it falls to about 2×10⁻¹¹, so the residual is entirely the
error-function approximation.

---

## Building and running on macOS (Tahoe / macOS 26)

### 1. Get gfortran

macOS doesn't ship a Fortran compiler, so you need gfortran first.
Two good options:

**Option A — Homebrew (recommended, easiest to keep updated)**

```bash
# Install Homebrew if you don't already have it:
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# Then install GCC, which includes gfortran:
brew install gcc
```

Homebrew will prompt you to install Apple's Command Line Tools first
if you don't have them — let it. After installation, `gfortran`
should be on your `PATH`; check with:

```bash
gfortran --version
```

If `gfortran` isn't found, Homebrew may have installed a versioned
binary like `gfortran-15` — either use that name directly, or add
`$(brew --prefix)/bin` to your `PATH` in `~/.zshrc`.

**Option B — Standalone installer**

fxcoudert publishes a standalone gfortran `.dmg` built specifically
for macOS Tahoe on Apple Silicon:
https://github.com/fxcoudert/gfortran-for-macOS/releases

(As of this writing there's no separate Intel build for Tahoe — if
you're on an Intel Mac, use Option A instead.) Download the
`.dmg`, run the installer, then confirm with `gfortran --version` in
a new terminal window.

### 2. Unzip/copy this folder to your Mac

Put all the files listed above in one folder, e.g.
`~/ibp_f77/`, and `cd` into it in Terminal.

### 3. Build

```bash
make
```

This compiles `pbutils.f`, `pbbatch.f`, and `pbsingle.f`, and links
two executables: `pbbatch` and `pbsingle`. If you'd rather not use
`make`, the equivalent manual commands are:

```bash
gfortran -std=legacy -O2 -c pbutils.f
gfortran -std=legacy -O2 -c pbbatch.f
gfortran -std=legacy -O2 -c pbsingle.f
gfortran -std=legacy -O2 -o pbbatch  pbutils.o pbbatch.o
gfortran -std=legacy -O2 -o pbsingle pbutils.o pbsingle.o
```

(`-std=legacy` tells gfortran to accept old-style F77 source without
warning about it; it's not required for correctness here, just to
keep the build quiet.)

### 4. Run

**Batch mode** (reads a file, writes a file):

```bash
./pbbatch
```
or
```bash
make run
```

The program reads `batch_inputs.txt` and `ibp_coefficients.txt` from
the current directory, and writes `loopresults.txt`. All three file
names are set near the top of `pbbatch.f` (`INFILE`, `OUTFILE`,
`CFFILE`) if you want to change them — F77 has no
command-line-argument handling worth relying on across compilers, so
this is the simplest robust way to point at different files.

**Interactive mode** (one run at a time, printed to the screen):

```bash
./pbsingle
```
or
```bash
make run-single
```

It will prompt you for DOY, longitude, F10.7, and local time in
turn, then print the full calculation — including all the
intermediate `DEBUG:` values and a coefficient sanity check — ending
with the occurrence probability. It only reads `ibp_coefficients.txt`;
no other files are touched.

### 5. Input file format (batch mode)

`batch_inputs.txt` is plain whitespace-separated numbers, one run per
line, with an **optional header row** (auto-detected and skipped):

```
DOY  Longitude(deg)  F10.7(sfu)  LocalTime(hr)
```

e.g. `80   30.0  120.0  20.50`. Add as many rows as you like.

### Cleaning up

```bash
make clean
```

removes the compiled `.o` files, the `pbbatch` and `pbsingle`
executables, and `loopresults.txt`.

### Troubleshooting: linker errors on Tahoe

If `gfortran` compiles the `.f` files fine (producing `.o` files)
but fails at the final link step with an error like

```
ld: unsupported tapi file type '!tapi-tbd' in YAML file ... libSystem.tbd
```

the cause is that the standalone gfortran build ships its own older
`ld`, which cannot parse the newer text-based library stubs
(`.tbd`) in the macOS SDK that Tahoe installs. Pointing `-isysroot`
at an older SDK does **not** help if that SDK also ships the new
`.tbd` format (it usually does on Tahoe).

The reliable fix is to compile with `gfortran` as usual but do the
**link step with the system `clang`**, which uses the current macOS
`ld` that understands the new stubs. Point it at gfortran's runtime
libraries and bake in an `-rpath` so the executable finds
`libgfortran` at run time:

```bash
GFLIB=$(dirname "$(gfortran -print-file-name=libgfortran.dylib)")
gfortran -std=legacy -O2 -c pbutils.f pbbatch.f pbsingle.f
clang -o pbbatch  pbutils.o pbbatch.o  -L"$GFLIB" -lgfortran -lquadmath -lm -Wl,-rpath,"$GFLIB"
clang -o pbsingle pbutils.o pbsingle.o -L"$GFLIB" -lgfortran -lquadmath -lm -Wl,-rpath,"$GFLIB"
```

This produces self-contained executables (no `DYLD_LIBRARY_PATH`
needed) and has been verified to reproduce the reference output
exactly on macOS Tahoe. You can ignore the harmless
`object file was built for newer 'macOS' version` link warnings.
