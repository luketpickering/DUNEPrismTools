#!/bin/bash

#Test variations

# for i in nu nubar; do
#   for j in X Y; do
#     for k in Horn1 Horn2 DecayPipe; do
#      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
#       -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${k}${j}Shift_${i}mode.mac \
#       -p test_alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}Shift_${i}/neutrino \
#       --number-of-jobs 1 \
#       --expected-disk 1GB \
#       --expected-mem 1GB \
#       --expected-walltime 4h
#     done
#   done
# done

for i in nu nubar; do
  for j in p1; do
    for k in WL HC DPR; do
     ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
      -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode_${k}${j}.mac \
      -p test_alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i} \
      --number-of-jobs 1000 \
      --expected-disk 1GB \
      --expected-mem 1GB \
      --expected-walltime 5h
    done
  done
done

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
#   -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_numode.mac \
#   -p nominal_1E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/neutrino \
#   --number-of-jobs 1000 --PPFX \
#   --expected-disk 1GB \
#   --expected-mem 1GB \
#   --expected-walltime 4h
#
${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
  -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_nubarmode.mac \
  -p nominal_2.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/antineutrino \
  --number-of-jobs 1500 --PPFX \
  --expected-disk 1GB \
  --expected-mem 1GB \
  --expected-walltime 5h
