#!/bin/bash

DO_PPFX="1"
PPFX_DIR="nominal_5E8POT_wppfx"
DO_PPFX_VARIATIONS="1"

DO_PPFX_COMPONENT_VARIATIONS="0"
PPFX_COMP_DIR="nominal_2.5E8POT_wallppfx"

DO_FOCUS="1"
FOCUS_DIR="Focussing"

DO_ALIGN="1"
ALIGN_DIR="Alignment"

for i in nu nubar; do

  if [ "${DO_PPFX}" == "1" ]; then
    if [ -e /pnfs/dune/persistent/users/picker24/${PPFX_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i} ]; then
      echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/picker24/${PPFX_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i}"
      continue
    fi

    if [ ! -e ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode.mac ]; then
      echo "[INFO]: No such g4lbnf macro OptimizedEngineeredNov2017Review_${i}mode.mac Not running."
      continue
    fi

    PPFX_ARG=""
    if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
      PPFX_ARG="--PPFX"
    fi

    ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
      -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode.mac \
      -p ${PPFX_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i} \
      --number-of-jobs 3500 ${PPFX_ARG} \
      --expected-disk 1GB \
      --expected-mem 1.9GB \
      --expected-walltime 4h
  fi

  if [ "${DO_PPFX_COMPONENT_VARIATIONS}" == "1" ]; then
    if [ -e /pnfs/dune/persistent/users/picker24/${PPFX_COMP_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i} ]; then
      echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/picker24/${PPFX_COMP_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i}"
      continue
    fi

    if [ ! -e ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode.mac ]; then
      echo "[INFO]: No such g4lbnf macro OptimizedEngineeredNov2017Review_${i}mode.mac Not running."
      continue
    fi

    ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
      -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode.mac \
      -p ${PPFX_COMP_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i} \
      --number-of-jobs 3500 --PPFX --Use-All-PPFX-Weights \
      --expected-disk 1GB \
      --expected-mem 1.9GB \
      --expected-walltime 4h
  fi

  #focussing
  if [ "${DO_FOCUS}" == "1" ]; then
    for j in p1 m1; do
      for k in WL HC DPR TargetDensity BeamSigma BeamOffsetX BeamTheta BeamThetaPhi; do

      if [ -e /pnfs/dune/persistent/users/picker24/${FOCUS_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i} ]; then
        echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/picker24/${FOCUS_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}"
        continue
      fi

      if [ ! -e ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode_${k}${j}.mac ]; then
        echo "[INFO]: No such g4lbnf macro OptimizedEngineeredNov2017Review_${i}mode_${k}${j}.mac Not running."
        continue
      fi

      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
       -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode_${k}${j}.mac \
       -p ${FOCUS_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i} \
       --number-of-jobs 2500 \
       --expected-disk 1GB \
       --expected-mem 1GB \
       --expected-walltime 4h
     done
   done
 fi

#alignment
  if [ "${DO_ALIGN}" == "1" ]; then
    for j in Horn1 Horn2; do
      for k in X Y XNeg X3mm XNeg3mm; do

      if [ -e /pnfs/dune/persistent/users/picker24/${ALIGN_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i} ]; then
        echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/picker24/${ALIGN_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i}"
        continue
      fi

      if [ ! -e ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${j}${k}Shift_${i}mode.mac ]; then
        echo "[INFO]: No such g4lbnf macro OptimizedEngineeredNov2017Review_${j}${k}Shift_${i}mode.mac Not running."
        continue
      fi

      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
       -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${j}${k}Shift_${i}mode.mac \
       -p ${ALIGN_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i} \
       --number-of-jobs 2500 \
       --expected-disk 1GB \
       --expected-mem 1GB \
       --expected-walltime 4h
      done
    done
  fi
done
