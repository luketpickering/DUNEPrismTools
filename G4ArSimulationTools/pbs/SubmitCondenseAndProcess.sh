#!/bin/bash

if ! mkdir subdir_FHC; then
  echo "[ERROR]: Failed to make submission dir \"subdir_FHC\", does it already exist?"
  exit 1
fi

cd subdir_FHC

${DUNEPRISMTOOLSROOT}/scripts/FarmCondenseAndProcess.sh \
  -G /mnt/scratch/calcuttj/DunePrism/FHC/ \
  -R /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker/FHC/ \
  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis \
  -f

cd -
