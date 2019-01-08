# !/bin/bash

dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/nu/dk2nulite/*.root" --fhicl build_FD_fluxes_ppfx_oldbin.fcl -o FD_nu_OptimizedEngineeredNov2017Review_old_uncert_binning_wppfx.root

dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/nubar/dk2nulite/*.root" --fhicl build_FD_fluxes_ppfx_oldbin.fcl -o FD_nubar_OptimizedEngineeredNov2017Review_old_uncert_binning_wppfx.root

for i in nu nubar; do
  for j in p1; do
    for k in WL HC DPR; do

    dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite/*.root" --fhicl build_FD_fluxes_oldbin.fcl -o FD_${i}_OptimizedEngineeredNov2017Review_old_uncert_binning_${k}${j}.root

    done
  done
done

for i in nu nubar; do
  for j in Horn1 Horn2; do
    for k in X Y; do
    dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/Alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i}/dk2nulite/*.root" --fhicl build_FD_fluxes_oldbin.fcl -o FD_${i}_OptimizedEngineeredNov2017Review_old_uncert_binning_${j}${k}Shift.root

    done
  done
done
