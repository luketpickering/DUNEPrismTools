#!/bin/bash

#if ! mkdir subdir_FHC; then
#  echo "[ERROR]: Failed to make submission dir \"subdir_FHC\", does it already exist?"
#  exit 1
#fi

#cd subdir_FHC

#${DUNEPRISMTOOLSROOT}/scripts/sim_scripts/g4py/FarmCondenseAndProcess.sh \
#  -G /mnt/scratch/calcuttj/DunePrism/FHC/G4Py_2018.3.13 \
#  -R /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker_withmec/FHC \
#  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC \
#  -f -P 5E16

#cd -

# if ! mkdir subdir_RHC; then
#   echo "[ERROR]: Failed to make submission dir \"subdir_RHC\", does it already exist?"
#   exit 1
# fi

# cd subdir_RHC

# ${DUNEPRISMTOOLSROOT}/scripts/sim_scripts/g4py/FarmCondenseAndProcess.sh \
#   -G /mnt/scratch/calcuttj/DunePrism/RHC/G4Py_2018.3.14 \
#   -R /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker_withmec/RHC \
#   -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/RHC \
#   -f -P 3.5E16

# cd -

if ! mkdir subdir_FHC_fd; then
  echo "[ERROR]: Failed to make submission dir \"subdir_FHC_fd\", does it already exist?"
  exit 1
fi

cd subdir_FHC_fd

${DUNEPRISMTOOLSROOT}/scripts/sim_scripts/g4py/FarmCondenseAndProcess.sh \
  -G /mnt/scratch/calcuttj/DunePrism/FHC_FD/G4Py_2018.3.19/ \
  -R /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker_withmec/FHC_FD \
  -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_FD \
  -C ${DUNEPRISMTOOLSROOT}/configs/run_plans/RunPlan.LArDUNEFD.25.8mx23.8mx56.6mFV.xml \
  -F "-nx 1 -dmn -1340,-1180,-2880 -dmx 1340,1180,2880 -V 50,50,50 -nt 5000 -T 1000" \
  -f -P 5E20

cd -

# if ! mkdir subdir_RHC_fd; then
#   echo "[ERROR]: Failed to make submission dir \"subdir_RHC_fd\", does it already exist?"
#   exit 1
# fi

# cd subdir_RHC_fd

# ${DUNEPRISMTOOLSROOT}/scripts/sim_scripts/g4py/FarmCondenseAndProcess.sh \
#   -G /mnt/scratch/calcuttj/DunePrism/RHC_FD/G4Py_2018.3.19/ \
#   -R /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker_withmec/RHC_FD \
#   -o /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/RHC_FD \
#   -C ${DUNEPRISMTOOLSROOT}/configs/run_plans/RunPlan.LArDUNEFD.25.8mx23.8mx56.6mFV.xml \
#   -F "-nx 1 -dmn -1340,-1180,-2880 -dmx 1340,1180,2880 -fv 50,50,50 -nt 5000 -T 1000" \
#   -f -P 5E20

# cd -


# if ! mkdir subdir_FHC2; then
#   echo "[ERROR]: Failed to make submission dir \"subdir_FHC2\", does it already exist?"
#   exit 1
# fi

# cd subdir_FHC2

# ${DUNEPRISMTOOLSROOT}/scripts/sim_scripts/g4py/FarmCondenseAndProcess.sh \
#   -G /mnt/scratch/calcuttj/DunePrism/FHC/G4Py_2018.3.18 \
#   -R /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker_withmec/FHC_2E8POT \
#   -o /mnt/home/picker24/Generation/GENIE/DUNE-PRISM/Analysis/FHC_2E8POT_p1 \
#   -f -P 5E16 -N 1000

# cd -

#if ! mkdir subdir_FHC3; then
#  echo "[ERROR]: Failed to make submission dir \"subdir_FHC3\", does it already exist?"
#  exit 1
#fi

#cd subdir_FHC3

#${DUNEPRISMTOOLSROOT}/scripts/FarmCondenseAndProcess.sh \
#  -G /mnt/scratch/calcuttj/DunePrism/FHC/G4Py_2018.3.18 \
#  -R /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker_withmec/FHC_2E8POT \
#  -o /mnt/home/picker24/Generation/GENIE/DUNE-PRISM/Analysis/FHC_2E8POT_p2 \
#  -f -P 5E16 -S 1000 -N 1000

#cd -
