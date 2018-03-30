#!/bin/bash

if ! mkdir subdir_FHCp; then
  echo "[ERROR]: Failed to make submission dir \"subdir_FHCp\", does it already exist?"
  exit 1
fi

cd subdir_FHCp

${DUNEPRISMTOOLSROOT}/scripts/FarmProcess.sh \
  -C /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Condensed.2018-03-20 \
  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps \
  -R ${DUNEPRISMTOOLSROOT}/configs/run_plans/RunPlan.39mLAr.4mx3mx5mActive_overlaps.xml \
  -P 5E16 -f

cd -

# if ! mkdir subdir_FHCp2; then
#   echo "[ERROR]: Failed to make submission dir \"subdir_FHCp\", does it already exist?"
#   exit 1
# fi

# cd subdir_FHCp2

# ${DUNEPRISMTOOLSROOT}/scripts/FarmProcess.sh \
#   -C /mnt/home/picker24/Generation/GENIE/DUNE-PRISM/Analysis/FHC_2E8POT_p1/Condensed.2018-03-21 \
#   -o /mnt/home/picker24/Generation/GENIE/DUNE-PRISM/Analysis/FHC_2E8POT_p1_withintermediates \
#   -R ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m_withintermediates.xml \
#   -P 5E16 -f

# cd -

# if ! mkdir subdir_FHCp3; then
#   echo "[ERROR]: Failed to make submission dir \"subdir_FHCp\", does it already exist?"
#   exit 1
# fi
#
# cd subdir_FHCp3
#
# ${DUNEPRISMTOOLSROOT}/scripts/FarmProcess.sh \
#   -C /mnt/home/picker24/Generation/GENIE/DUNE-PRISM/Analysis/FHC_2E8POT_p2/Condensed.2018-03-21 \
#   -o /mnt/home/picker24/Generation/GENIE/DUNE-PRISM/Analysis/FHC_2E8POT_p2_withintermediates \
#   -R ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m_withintermediates.xml \
#   -P 5E16 -f
#
# cd -
