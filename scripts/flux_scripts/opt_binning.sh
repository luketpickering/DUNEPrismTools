#!/bin/bash
NFILES=100
NEVEQuiv_RS=50000
NEVEQuiv_WS=25000
NEVEQuiv_NUE=10000
NUSTUB=/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/nu/dk2nulite
NUBARSTUB=/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/nubar/dk2nulite

ls ${NUSTUB}/g4*.root | head -${NFILES} | sed "s|/pnfs/dune|root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune|g" > nu.files.list
ls ${NUBARSTUB}/g4*.root | head -${NFILES} | sed "s|/pnfs/dune|root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune|g" > nubar.files.list

#cat nu.files.list
#cat nubar.files.list

#exit 0

for hc in nu nubar; do
  for i in 14,numu -14,numubar 12,nue -12,nuebar; do
    PDG=$(echo ${i} | cut -d , -f 1)
    NAME=$(echo ${i} | cut -d , -f 2)

    if [ ${hc} == "nu" ] && [ ${PDG} == 14 ]; then
      NEVEQuiv=${NEVEQuiv_RS}
    elif [ ${hc} == "nu" ] && [ ${PDG} == -14 ]; then
      NEVEQuiv=${NEVEQuiv_WS}
    elif [ ${hc} == "nubar" ] && [ ${PDG} == 14 ]; then
      NEVEQuiv=${NEVEQuiv_WS}
    elif [ ${hc} == "nubar" ] && [ ${PDG} == -14 ]; then
      NEVEQuiv=${NEVEQuiv_RS}
    elif [ ${PDG} == 12 ] || [ ${PDG} == -12 ]; then
      NEVEQuiv=${NEVEQuiv_NUE}
    fi

    if ! dp_BuildFluxes \
          --fhicl build_ND_fluxes_${hc}mode_dummybin_offaxis.fcl \
          -i ${hc}.files.list \
          --only-pdg ${PDG} \
          --optimized-binning-file ${hc}.${NAME}.optbin.fcl \
          --optimized-binning-N ${NEVEQuiv} \
          --optimized-binning-min-width 0.025 \
          --optimized-binning-max-width 15.0; then 
      exit 1 
    fi

  done
done
