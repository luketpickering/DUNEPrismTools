#!/bin/bash

if ! mkdir subdir_FHCp; then
  echo "[ERROR]: Failed to make submission dir \"subdir_FHCp\", does it already exist?"
  exit 1
fi

cd subdir_FHCp

${DUNEPRISMTOOLSROOT}/scripts/FarmProcess.sh \
  -C /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis_withmec/FHC/Condensed.2018-03-13 \
  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis_withmec/FHC \
  -R ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -P 5E16 -f

cd -
