# Makefile for the IBP single-run model sharing package.
# Tested with GNU Fortran (gfortran) on macOS and Linux.
#
# Two ways to use this model:
#   1. pbsingle    - standalone interactive program (prompts for
#                    input, prints full debug detail). Good for
#                    manual testing / demos.
#   2. libibp.a    - a linkable static library exposing IBP_INIT and
#                    IBP_PROBABILITY, with no I/O and no STOP calls,
#                    for embedding directly in another model.
#                    test_ibp_probability.f shows how to call it.

FC      = gfortran
FFLAGS  = -std=legacy -O2

all: pbsingle libibp.a test_ibp

# --- Standalone interactive program ---
pbsingle: pbsingle.o pbutils.o
	$(FC) $(FFLAGS) -o pbsingle pbsingle.o pbutils.o

pbsingle.o: pbsingle.f pbcoeff.inc
	$(FC) $(FFLAGS) -c pbsingle.f

# --- Callable library ---
libibp.a: ibp_probability.o pbutils.o
	ar rcs libibp.a ibp_probability.o pbutils.o

ibp_probability.o: ibp_probability.f pbcoeff.inc
	$(FC) $(FFLAGS) -c ibp_probability.f

pbutils.o: pbutils.f pbcoeff.inc
	$(FC) $(FFLAGS) -c pbutils.f

# --- Test harness / usage example for the library ---
test_ibp: test_ibp_probability.o libibp.a
	$(FC) $(FFLAGS) -o test_ibp test_ibp_probability.o -L. -libp

test_ibp_probability.o: test_ibp_probability.f
	$(FC) $(FFLAGS) -c test_ibp_probability.f

test: test_ibp
	./test_ibp

run: pbsingle
	./pbsingle

clean:
	rm -f *.o *.a pbsingle test_ibp

.PHONY: all run test clean
