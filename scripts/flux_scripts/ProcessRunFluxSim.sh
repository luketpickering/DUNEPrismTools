#!/bin/bash

for i in nu nubar; do

    if [ -e /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i} ]; then
      echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i}"
      continue
    fi

    if [ ! -e ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode.mac ]; then
      echo "[INFO]: No such g4lbnf macro OptimizedEngineeredNov2017Review_${i}mode.mac Not running."
      continue
    fi

    ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
      -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode.mac \
      -p nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i} \
      --number-of-jobs 3500 --PPFX \
      --expected-disk 1GB \
      --expected-mem 1.9GB \
      --expected-walltime 4h

done

#focussing
for i in nu nubar; do
 for j in p1 m1; do
   for k in WL HC DPR TargetDensity BeamSigma BeamOffsetX BeamTheta BeamThetaPhi; do

    if [ -e /pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i} ]; then
      echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}"
      continue
    fi

    if [ ! -e ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode_${k}${j}.mac ]; then
      echo "[INFO]: No such g4lbnf macro OptimizedEngineeredNov2017Review_${i}mode_${k}${j}.mac Not running."
      continue
    fi

    ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
     -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode_${k}${j}.mac \
     -p Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i} \
     --number-of-jobs 2500 \
     --expected-disk 1GB \
     --expected-mem 1GB \
     --expected-walltime 4h
   done
 done
done

#alignment
for i in nu nubar; do
  for j in Horn1 Horn2; do
    for k in X Y XNeg X3mm XNeg3mm; do

    if [ -e /pnfs/dune/persistent/users/picker24/Alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i} ]; then
      echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/picker24/Alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i}"
      continue
    fi

    if [ ! -e ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${j}${k}Shift_${i}mode.mac ]; then
      echo "[INFO]: No such g4lbnf macro OptimizedEngineeredNov2017Review_${j}${k}Shift_${i}mode.mac Not running."
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
