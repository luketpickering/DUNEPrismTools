# !/bin/bash

dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/nominal_1.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/neutrino/dk2nulite/*.root" --fhicl build_FD_fluxes_numode_ppfx.fcl -o FD_nu_OptimizedEngineeredNov2017Review_uncert_binning_wppfx.root

dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/antineutrino/dk2nulite/*.root" --fhicl build_FD_fluxes_nubarmode_ppfx.fcl -o FD_nubar_OptimizedEngineeredNov2017Review_uncert_binning_wppfx.root

for i in nu nubar; do
  for j in p1; do
    for k in WL HC DPR; do

    dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite/*.root" --fhicl build_FD_fluxes_${i}mode_ppfx.fcl -o FD_${i}_OptimizedEngineeredNov2017Review_uncert_binning_${k}${j}.root

    done
  done
done
