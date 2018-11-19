#!/bin/bash

dp_CombineBuiltFluxes  \
	-i "/pnfs/dune/persistent/users/picker24/nominal_1.5E8_wppfx/DUNEPrismFluxes/ND_nu/uncert_binning/flux/Fluxes.*.root" \
  --NPPFXU 100 \
	-o ND_nu_OptimizedEngineeredNov2017Review_uncert_binning_wppfx.root

dp_CombineBuiltFluxes  \
	-i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8_wppfx/DUNEPrismFluxes/ND_nubar/uncert_binning/flux/Fluxes.*.root" \
  --NPPFXU 100 \
	-o ND_nubar_OptimizedEngineeredNov2017Review_uncert_binning_wppfx.root

#With focussing
for i in nu nubar; do
  for j in p1; do
    for k in WL HC DPR; do

		dp_CombineBuiltFluxes  \
			-i "/pnfs/dune/persistent/users/picker24/Focussing/DUNEPrismFluxes/ND_${i}/${k}${j}/uncert_binning/flux/Fluxes.*.root" \
			-o ND_${i}_OptimizedEngineeredNov2017Review_uncert_binning_${k}${j}.root

		done
	done
done
