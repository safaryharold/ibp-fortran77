# Makefile for the FORTRAN 77 IBP (Ionospheric Bubble Probability)
# model. Tested with GNU Fortran (gfortran) on macOS and Linux.
#
# Two programs share the same coefficient loader / physics
# functions (pbutils.f) and COMMON block (pbcoeff.inc):
#   pbbatch   - batch driver, reads batch_inputs.txt, writes
#               loopresults.txt
#   pbsingle  - interactive, single-run driver; prompts for one
#               set of inputs and prints the result with debug
#               detail

FC      = gfortran
FFLAGS  = -std=legacy -O2
LDFLAGS =

all: pbbatch pbsingle

pbbatch: pbbatch.o pbutils.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o pbbatch pbbatch.o pbutils.o

pbsingle: pbsingle.o pbutils.o
	$(FC) $(FFLAGS) $(LDFLAGS) -o pbsingle pbsingle.o pbutils.o

pbbatch.o: pbbatch.f pbcoeff.inc
	$(FC) $(FFLAGS) -c pbbatch.f

pbsingle.o: pbsingle.f pbcoeff.inc
	$(FC) $(FFLAGS) -c pbsingle.f

pbutils.o: pbutils.f pbcoeff.inc
	$(FC) $(FFLAGS) -c pbutils.f

run: pbbatch
	./pbbatch

run-single: pbsingle
	./pbsingle

clean:
	rm -f *.o pbbatch pbsingle loopresults.txt

.PHONY: all run run-single clean
