#!/bin/bash

# #focussing
# for i in nu nubar; do
#  for j in p1; do
#    for k in WL HC DPR; do
#     ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
#      -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode_${k}${j}.mac \
#      -p Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i} \
#      --number-of-jobs 2500 \
#      --expected-disk 1GB \
#      --expected-mem 1GB \
#      --expected-walltime 4h
#    done
#  done
# done

#alignment
for i in nu nubar; do
  for j in Horn1 Horn2; do
    for k in X Y; do
# for i in nu; do
#   for j in Horn1; do
#     for k in X; do

    if [ -e /pnfs/dune/persistent/users/picker24/Alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i} ]; then
      echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/picker24/Alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i}"
      continue
    fi

    ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
     -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${j}${k}Shift_${i}mode.mac \
     -p Alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i} \
     --number-of-jobs 2500 \
     --expected-disk 1GB \
     --expected-mem 1GB \
     --expected-walltime 4h
   done
 done
done

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
#   -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_numode.mac \
#   -p nominal_3.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/neutrino \
#   --number-of-jobs 3500 --PPFX \
#   --expected-disk 1GB \
#   --expected-mem 1.9GB \
#   --expected-walltime 4h
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
#   -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_nubarmode.mac \
#   -p nominal_3.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/antineutrino \
#   --number-of-jobs 3500 --PPFX \
#   --expected-disk 1GB \
#   --expected-mem 1.9GB \
#   --expected-walltime 4h
