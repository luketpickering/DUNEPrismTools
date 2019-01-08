#!/bin/bash

# dp_CombineBuiltFluxes  \
# 	-i "/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/DUNEPrismFluxes/ND_nu/uncert_binning/flux/Fluxes.*.root" \
#   --NPPFXU 100 \
# 	-o ND_nu_OptimizedEngineeredNov2017Review_uncert_binning_wppfx.root
#
# dp_CombineBuiltFluxes  \
# 	-i "/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/DUNEPrismFluxes/ND_nubar/uncert_binning/flux/Fluxes.*.root" \
#   --NPPFXU 100 \
# 	-o ND_nubar_OptimizedEngineeredNov2017Review_uncert_binning_wppfx.root

#With focussing
# for i in nu nubar; do
#   for j in p1; do
#     for k in WL HC DPR; do
#
# 		dp_CombineBuiltFluxes  \
# 			-i "/pnfs/dune/persistent/users/picker24/Focussing/DUNEPrismFluxes/ND_${i}/${k}${j}/uncert_binning/flux/Fluxes.*.root" \
# 			-o ND_${i}_OptimizedEngineeredNov2017Review_uncert_binning_${k}${j}.root
#
# 		done
# 	done
# done

# for i in nu nubar; do
#   for j in DecayPipe Horn1 Horn2; do
#     for k in X Y; do
for i in nu ; do
  for j in Horn1; do
    for k in X; do


			dp_CombineBuiltFluxes  \
				-i "/pnfs/dune/persistent/users/picker24/Alignment/DUNEPrismFluxes/ND_${i}/${j}${k}/uncert_binning/flux/Fluxes.*.root" \
				-o ND_${i}_OptimizedEngineeredNov2017Review_uncert_binning_${j}${k}Shift.root

    done
  done
done
