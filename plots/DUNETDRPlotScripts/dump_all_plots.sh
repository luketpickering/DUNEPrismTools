#!/bin/bash

if ! hash root-config; then
  echo "[ERROR]: ROOT is not set up."
  exit 1
fi

if [ ! -e TH2Jagged/build/Linux/lib/libTH2Jagged.so ]; then
  rm -rf TH2Jagged
  git clone https://github.com/luketpickering/TH2Jagged.git
  cd TH2Jagged; mkdir build; cd build; cmake ../; make install
  cd ../../;
fi

rm -rf tikz tikz_render plots

mkdir tikz tikz_render

root -l -b -q ComponentFluxPreds.C
root -l -b -q PlotNearFarRatioErrorSources.C
root -l -b -q PlotNearOffAxisErrorSources.C
root -l -b -q DumpBinningPlots.C
root -l -b -q DumpMatrix.C

./render.tikz.sh tikz tikz_render

mkdir plots
mv *.pdf tikz_render/*.pdf plots/
