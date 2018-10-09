#!/bin/bash

# if ! mkdir subdir_FHCp; then
#   echo "[ERROR]: Failed to make submission dir \"subdir_FHCp\", does it already exist?"
#   exit 1
# fi
#
# cd subdir_FHCp
#
# ${DUNEPRISMTOOLSROOT}/scripts/sim_scripts/g4py/FarmProcess.sh \
#   -C /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Condensed.2018-03-30 \
#   -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps \
#   -R ${DUNEPRISMTOOLSROOT}/configs/run_plans/RunPlan.39mLAr.4mx3mx5mActive_overlaps.xml \
#   -P 5E16 -f
#
# cd -

if ! mkdir subdir_FHCp7m; then
  echo "[ERROR]: Failed to make submission dir \"subdir_FHCp7m\", does it already exist?"
  exit 1
fi

cd subdir_FHCp7m

${DUNEPRISMTOOLSROOT}/scripts/sim_scripts/g4py/FarmProcess.sh \
  -C /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Condensed.2018-03-30 \
  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m \
  -R ${DUNEPRISMTOOLSROOT}/configs/run_plans/RunPlan.39mLAr.7mx3mx5mActive_overlaps.xml \
  -P 5E16 -f

cd -
