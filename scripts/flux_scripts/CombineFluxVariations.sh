#!/bin/bash

DET="FD"
IDIR="old_uncert_binning"

for i in nu nubar; do

	dp_CombineBuiltFluxes  \
		-i "/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/DUNEPrismFluxes/${DET}_${i}/${IDIR}/flux/Fluxes.*.root" \
	  --NPPFXU 100 \
		-o ${DET}_${i}_OptimizedEngineeredNov2017Review_${IDIR}_wppfx.root

done

#With focussing
for i in nu nubar; do
  for j in p1 m1; do
    for k in WL HC DPR; do

      if [ ! -e /pnfs/dune/persistent/users/picker24/Focussing/DUNEPrismFluxes/${DET}_${i}/${k}${j}/${IDIR}/flux ]; then
        continue;
      fi

		dp_CombineBuiltFluxes  \
			-i "/pnfs/dune/persistent/users/picker24/Focussing/DUNEPrismFluxes/${DET}_${i}/${k}${j}/${IDIR}/flux/Fluxes.*.root" \
			-o ${DET}_${i}_OptimizedEngineeredNov2017Review_${IDIR}_${k}${j}.root

		done
	done
done

for i in nu nubar; do
  for j in Horn1 Horn2; do
    for k in X Y XNeg; do

      if [ ! -e /pnfs/dune/persistent/users/picker24/Alignment/DUNEPrismFluxes/${DET}_${i}/${j}${k}/${IDIR}/flux ]; then
        continue;
      fi

			dp_CombineBuiltFluxes  \
				-i "/pnfs/dune/persistent/users/picker24/Alignment/DUNEPrismFluxes/${DET}_${i}/${j}${k}/${IDIR}/flux/Fluxes.*.root" \
				-o ${DET}_${i}_OptimizedEngineeredNov2017Review_${IDIR}_${j}${k}Shift.root

    done
  done
done
