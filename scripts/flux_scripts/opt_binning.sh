#!/bin/bash
NFILES=20
NEVEQuiv=100000
NUSTUB=/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/nubar/dk2nulite
NUBARSTUB=/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/nubar/dk2nulite

ls ${NUSTUB}/g4*.root | head -${NFILES} | sed "s|/pnfs/dune|root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune|g" > nu.files.list
ls ${NUBARSTUB}/g4*.root | head -${NFILES} | sed "s|/pnfs/dune|root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune|g" > nubar.files.list

for hc in nu nubar; do
  for i in 14,numu -14,numubar 12,nue -12,nuebar; do
    PDG=$(echo ${i} | cut -d , -f 1)
    NAME=$(echo ${i} | cut -d , -f 1)

    dp_BuildFluxes --fhicl build_ND_fluxes_${hc}mode_dummybin_offaxis.fcl -i ${hc}.files.list --only-pdg ${PDG} --optimized-binning-file ${hc}.${NAME}.optbin.fcl --optimized-binning-N ${NEVEQuiv} --optimized-binning-min-width 0.025

  done
done
