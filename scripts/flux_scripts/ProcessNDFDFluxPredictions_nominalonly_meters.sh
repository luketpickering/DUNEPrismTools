# !/bin/bash

SUBDIR="fine_binning_meters"
FORCEOVERWRITE="true"

EDISK="2GB"
ETIME_PPFX="2h"
ETIME="1h"

for i in nu nubar; do

    if [ ! -e /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i}/dk2nulite ]; then
      echo "[INFO]: No input directory for ${i} wppfx, skipping."
      continue
    fi

    if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/DUNEPrismFluxes/ND_${i}/${SUBDIR} ]; then
      echo "[INFO]: Already have ${i} wppfx not reprocessing."
      continue
    fi

  ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
     --expected-walltime ${ETIME_PPFX} --expected-disk ${EDISK} \
     --expected-mem 512MB -FN build_ND_fluxes_fitbinning_meters.fcl \
     -FF build_FD_fluxes_fitbinning_meters.fcl \
     -p nominal_5E8POT_wppfx/DUNEPrismFluxes/__DET___${i}/${SUBDIR} \
     -i /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i}/dk2nulite \
     -n 20 -f

done
