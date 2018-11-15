# !/bin/bash

BINNING_DESCRIPTOR_FITS="0_10:0.025"
BINNING_DESCRIPTOR_UNCERTS="0,0.5,1_3:0.25,3_4:0.5,4_10:1,10_20:2"

dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/nominal_1.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/neutrino/dk2nulite/*.root"
-x -12.90_12.90:25.80 -h 2260 -vb ${BINNING_DESCRIPTOR_FITS} -z 128700000 -L --PPFX -o /pnfs/dune/persistent/users/picker24/nominal_1.5E8_wppfx/DUNEPrismFluxes/FD_nu/fit_binning/FD_nu_OptimizedEngineeredNov2017Review_fit_binning_wppfx.root

dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/nominal_1.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/neutrino/dk2nulite/*.root"
-x -12.90_12.90:25.80 -h 2260 -vb ${BINNING_DESCRIPTOR_UNCERTS} -z 128700000 -L --PPFX -o /pnfs/dune/persistent/users/picker24/nominal_1.5E8_wppfx/DUNEPrismFluxes/FD_nu/uncert_binning/FD_nu_OptimizedEngineeredNov2017Review_uncert_binning_wppfx.root

dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/antineutrino/dk2nulite/*.root"
-x -12.90_12.90:25.80 -h 2260 -vb ${BINNING_DESCRIPTOR_FITS} -z 128700000 -L --PPFX -o /pnfs/dune/persistent/users/picker24/nominal_2.5E8_wppfx/DUNEPrismFluxes/FD_nubar/fit_binning/FD_nubar_OptimizedEngineeredNov2017Review_fit_binning_wppfx.root

dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/antineutrino/dk2nulite/*.root"
-x -12.90_12.90:25.80 -h 2260 -vb ${BINNING_DESCRIPTOR_UNCERTS} -z 128700000 -L --PPFX -o /pnfs/dune/persistent/users/picker24/nominal_2.5E8_wppfx/DUNEPrismFluxes/FD_nubar/uncert_binning/FD_nubar_OptimizedEngineeredNov2017Review_uncert_binning_wppfx.root


for i in nu nubar; do
  for j in p1; do
    for k in WL HC DPR; do

    dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite/*.root"
    -x -12.90_12.90:25.80 -h 2260 -vb ${BINNING_DESCRIPTOR_FITS} -z 128700000 -L --PPFX -o /pnfs/dune/persistent/users/picker24/Focussing/DUNEPrismFluxes/FD_${i}/${k}${j}/fit_binning/FD_${i}_OptimizedEngineeredNov2017Review_fit_binning_wppfx.root

    dp_BuildFluxes -i "/pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite/*.root"
    -x -12.90_12.90:25.80 -h 2260 -vb ${BINNING_DESCRIPTOR_UNCERTS} -z 128700000 -L --PPFX -o /pnfs/dune/persistent/users/picker24/Focussing/DUNEPrismFluxes/FD_${i}/${k}${j}/uncert_binning/FD_${i}_OptimizedEngineeredNov2017Review_uncert_binning_wppfx.root


    done
  done
done
