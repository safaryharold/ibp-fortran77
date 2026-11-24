#!/bin/sh
# Download, build and run IRI-2026 with the IBP plasma-bubble output.
#
# Usage:  ./run_iri2026_ibp.sh            (macOS + Anaconda: FC="gfortran -B/usr/bin/" ./run_iri2026_ibp.sh)
#
# Everything is placed in ./IRI-2026-build (ignored by git). The run uses
# iritest_ibp_input.txt: geographic lat 0, lon -76 (Jicamarca sector),
# 21 March 2024, 400 km, local time 0-23 h, output option 3 with
# parameter 92 = plasma bubble probability (OARR(92), jf(38)=.true.).
set -e
FC=${FC:-gfortran}
HERE=$(cd "$(dirname "$0")" && pwd)
B="$HERE/IRI-2026-build"
mkdir -p "$B"
cd "$B"

if [ ! -f irisub.for ]; then
  echo "Downloading IRI-2026 source from https://irimodel.org/IRI-2026/ ..."
  curl -fsSL -A "Mozilla/5.0" -O https://irimodel.org/IRI-2026/00_iri.tar
  tar xf 00_iri.tar
  mv IRI-zip/* . && rmdir IRI-zip
fi
for f in apf107.dat ig_rz.dat; do
  [ -f "$f" ] || curl -fsSL -A "Mozilla/5.0" -o "$f" "https://irimodel.org/indices/$f"
done

$FC -std=legacy -w -O1 -o iritest iritest.for irisub.for irifun.for \
    iritec.for iridreg.for igrf.for cira.for iriflip.for rocdrift.for

rm -f fort.7
./iritest < "$HERE/iritest_ibp_input.txt" > iritest_screen.txt 2>&1
echo "----- IRI-2026 output (fort.7) -----"
cat fort.7
