#!/bin/bash

if ! mkdir subdir_FHCp; then
  echo "[ERROR]: Failed to make submission dir \"subdir_FHCp\", does it already exist?"
  exit 1
fi

cd subdir_FHCp

${DUNEPRISMTOOLSROOT}/scripts/FarmProcess.sh \
  -C /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/Condensed.2018-03-05 \
  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis \
  -N 2

cd -
