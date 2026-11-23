Validation and reproducibility materials for the FORTRAN 77 IBP port.

DATA
  batch_inputs_504.txt         504-case input grid: DOY, longitude (deg E),
                               F10.7 (sfu), local time (hr).
  fortran_loopresults_504.txt  pbbatch output for that grid (after the
                               calendar and overflow-clamp fixes).

REFERENCE IMPLEMENTATION AND CROSS-VALIDATION
  ibp_reference.py   Independent, dependency-light Python implementation of
                     the same published model (Stolle et al., 2024; official
                     IBP release v2.0.2). Reads the same ibp_coefficients.txt
                     and uses the same closed-form (error-function)
                     probability expression, with exact math.erf.
  compare_504.py     Evaluates ibp_reference.py over batch_inputs_504.txt and
                     differences it row-by-row against
                     fortran_loopresults_504.txt. Reproduces the paper's
                     cross-validation:
                        max  |P_ref - P_fortran|  = 8.2e-7
                        mean |P_ref - P_fortran|  = 4.4e-8
                     The residual is dominated by the Fortran port's ERF77
                     error-function approximation (the reference uses exact
                     erf).  Run:  python compare_504.py

  The reference used in the paper is the official IBP v2.0.2 Python package
  (ibpmodel / ibpcalc.py) from https://igit.iap-kborn.de/ibp/ibp-model .
  ibp_reference.py is an independent re-implementation included so this
  comparison is reproducible without installing that package; it reproduces
  both the official model and the Fortran port to the ERF77 level above.

FIGURE REGENERATION
  make_figures.py    Regenerates the climatology figures directly from
                     ibp_reference.py:
                        figures/figure3_seasonal.png  (4 seasons, F10.7=150)
                        figures/figure4_sectors.png   (4 sectors, Dec, F10.7=150)
                        figures/figure5_doy_lt.png     (70 W, F10.7=150)
                        figures/figure6_monthly.png    (monthly mean, F10.7=80..200)
                        figures/figure7_crossval.png   (Python vs Fortran + difference)
                     Each figure is also written as a vector .pdf (embedded fonts)
                     next to the .png.
                     figure7 runs ../ibp_model_f77/pbbatch for the Fortran panel if
                     that binary has been built (otherwise it falls back to the Python
                     reference and the difference panel is zero by construction).
                     Requires numpy + matplotlib (seaborn optional, used for styling).
                     Run:  python make_figures.py
  make_figure1.py    Draws the model data-flow diagram (Figure 1) as vector
                     figures/Figure_1.pdf and .eps plus a 1000 dpi .png.
                     Requires matplotlib only.  Run:  python make_figure1.py

COMPARISON WITH THE OFFICIAL PYTHON PACKAGE
  ibp_dump.f           Tiny driver: reads "DOY LON F107 LT" lines, prints the
                       full-precision IBP_PROBABILITY result.
  compare_official.py  50,000 random inputs through ibp_dump and through the
                       official ibpmodel code, with the exact erf and with ERF77
                       substituted. Expected (v2.0.2 coefficients):
                         exact erf : max ~2.2e-06, mean ~3e-08
                         ERF77     : max ~2e-11
                       Needs: pip install ibpmodel==2.0.2, numpy, scipy, cdflib.
                       (ibpmodel 2.0 imports PyTorch for its optional
                       neural-network mode, which this test does not use.)
                       To test v1.4.0 instead: pip install ibpmodel==1.4.2,
                       convert its ..._0004.cdf with cdf_to_coefficients.py
                       into ibp_coefficients.txt (expected max ~2.0e-06).
                       Build:  (cd ../ibp_single_model_share && make libibp.a)
                               gfortran -std=legacy -O2 ibp_dump.f \
                                 -L../ibp_single_model_share -libp -o ibp_dump
                       Run:    python compare_official.py
  ibp_components.f     Driver that prints the Fortran values of g_osc, ts_month,
                       phi_month, lambda, mu, P and ERF77 for each input line.
  make_figure8.py      Figure 8: component-by-component and full-domain
                       comparison with the official ibpmodel 2.0.2 code over
                       the same 50,000 inputs as compare_official.py.
                       Build:  gfortran -std=legacy -O2 ibp_components.f \
                                 -I../ibp_single_model_share \
                                 -L../ibp_single_model_share -libp -o ibp_components
                       Run:    python make_figure8.py
                       Expected: model components agree to <=1.3e-15 (relative),
                       phi_month identical, ERF77 1.4e-07, P 2.2e-06.
  cdf_to_coefficients.py
                       Converts an official coefficient CDF (shipped inside the
                       ibpmodel package, e.g. ..._0101.cdf in v2.0.2) into the
                       ibp_coefficients.txt format read by the Fortran code.
                       Use it to run a newer IBP release without recompiling,
                       then pass the same CDF to compare_official.py.

REQUIREMENTS
  compare_504.py  : Python 3 standard library only.
  make_figures.py : Python 3 + numpy + matplotlib.
  make_figure1.py : Python 3 + matplotlib.
  ibp_coefficients.txt must sit alongside these scripts (a copy is included).
