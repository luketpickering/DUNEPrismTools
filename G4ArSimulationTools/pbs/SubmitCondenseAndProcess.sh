#!/bin/bash

# if ! mkdir subdir_FHC_fd; then
#   echo "[ERROR]: Failed to make submission dir \"subdir_FHC_fd\", does it already exist?"
#   exit 1
# fi

# cd subdir_FHC_fd

# FHC ND -P 5E16

# ${DUNEPRISMTOOLSROOT}/scripts/FarmCondenseAndProcess.sh \
#   -G /mnt/scratch/calcuttj/DunePrism/FHC_FD/G4Py_2018.3.13 \
#   -R /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker_withmec/FHC_FD \
#   -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis_withmec/FHC_FD \
#   -C ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml \
#   -F "-nx 3 -dmn -1340,-1180,-2880 -dmx 1340,1180,2880 -fv 50,50,50 -nt 5000 -T 1000" \
#   -f -N 2 -P 5E20

# cd -


if ! mkdir subdir_FHC_fd; then
  echo "[ERROR]: Failed to make submission dir \"subdir_FHC_fd\", does it already exist?"
  exit 1
fi

cd subdir_FHC_fd

${DUNEPRISMTOOLSROOT}/scripts/FarmCondenseAndProcess.sh \
  -G /mnt/scratch/calcuttj/DunePrism/FHC_FD/G4Py_2018.3.19/ \
  -R /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker_withmec/FHC_FD \
  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis_withmec/FHC_FD \
  -C ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml \
  -F "-nx 3 -dmn -1340,-1180,-2880 -dmx 1340,1180,2880 -fv 50,50,50 -nt 1000 -T 250" \
  -f -P 5E20

cd -

if ! mkdir subdir_RHC_fd; then
  echo "[ERROR]: Failed to make submission dir \"subdir_RHC_fd\", does it already exist?"
  exit 1
fi

cd subdir_RHC_fd

${DUNEPRISMTOOLSROOT}/scripts/FarmCondenseAndProcess.sh \
  -G /mnt/scratch/calcuttj/DunePrism/RHC_FD/G4Py_2018.3.19/ \
  -R /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker_withmec/RHC_FD \
  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis_withmec/RHC_FD \
  -C ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml \
  -F "-nx 3 -dmn -1340,-1180,-2880 -dmx 1340,1180,2880 -fv 50,50,50 -nt 1000 -T 250" \
  -f -P 5E20

cd -

#
# if ! mkdir subdir_RHC; then
#   echo "[ERROR]: Failed to make submission dir \"subdir_RHC\", does it already exist?"
#   exit 1
# fi

# cd subdir_RHC

# ${DUNEPRISMTOOLSROOT}/scripts/FarmCondenseAndProcess.sh \
#   -G /mnt/scratch/calcuttj/DunePrism/RHC/G4Py_2018.3.14 \
#   -R /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker_withmec/RHC \
#   -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis_withmec/RHC \
#   -f -P 3.5E16

# cd -
