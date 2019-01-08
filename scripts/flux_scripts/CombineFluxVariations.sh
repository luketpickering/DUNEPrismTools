#!/bin/bash

IDIR="old_uncert_binning"

dp_CombineBuiltFluxes  \
	-i "/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/DUNEPrismFluxes/ND_nu/${IDIR}/flux/Fluxes.*.root" \
  --NPPFXU 100 \
	-o ND_nu_OptimizedEngineeredNov2017Review_${IDIR}_wppfx.root

dp_CombineBuiltFluxes  \
	-i "/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/DUNEPrismFluxes/ND_nubar/${IDIR}/flux/Fluxes.*.root" \
  --NPPFXU 100 \
	-o ND_nubar_OptimizedEngineeredNov2017Review_${IDIR}_wppfx.root

#With focussing
for i in nu nubar; do
  for j in p1 m1; do
    for k in WL HC DPR; do

      if [ ! -e /pnfs/dune/persistent/users/picker24/Focussing/DUNEPrismFluxes/ND_${i}/${k}${j}/${IDIR}/flux ]; then
        continue;
      fi

		dp_CombineBuiltFluxes  \
			-i "/pnfs/dune/persistent/users/picker24/Focussing/DUNEPrismFluxes/ND_${i}/${k}${j}/${IDIR}/flux/Fluxes.*.root" \
			-o ND_${i}_OptimizedEngineeredNov2017Review_${IDIR}_${k}${j}.root

		done
	done
done

for i in nu nubar; do
  for j in Horn1 Horn2; do
    for k in X Y XNeg; do

      if [ ! -e /pnfs/dune/persistent/users/picker24/Focussing/DUNEPrismFluxes/ND_${i}/${k}${j}/${IDIR}/flux ]; then
        continue;
      fi

			dp_CombineBuiltFluxes  \
				-i "/pnfs/dune/persistent/users/picker24/Alignment/DUNEPrismFluxes/ND_${i}/${j}${k}/${IDIR}/flux/Fluxes.*.root" \
				-o ND_${i}_OptimizedEngineeredNov2017Review_${IDIR}_${j}${k}Shift.root

    done
  done
done
