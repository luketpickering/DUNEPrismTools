#!/bin/bash

IDIR="fine_binning"

for DET in "ND" "FD"; do
	for i in nu nubar; do
	dp_CombineBuiltFluxes  \
		-i "/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/DUNEPrismFluxes/${DET}_${i}/${IDIR}/flux/Fluxes.*.root" \
		-o ${DET}_${i}_OptimizedEngineeredNov2017Review_${IDIR}.root
	done
done
