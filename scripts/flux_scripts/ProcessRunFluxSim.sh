#!/bin/bash

if [ -z ${DUNEPRISMTOOLSROOT} ]; then
  echo "[ERROR]: DUNE-PRISM Tools environment is not set up."
  return;
fi

DO_PPFX="0"
PPFX_OUTPUT_DIR="nominal_5E8POT_wppfx"
DO_PPFX_VARIATIONS="1"

DO_HIGHERHC="1"
HIGHERHC_DIR="HigherHC_2.5E8POT_2019117"
if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
  HIGHERHC_DIR="HigherHC_2.5E8POT_wppfx_2019117"
fi

DO_PPFX_COMPONENT_VARIATIONS="0"
PPFX_COMP_DIR="nominal_2.5E8POT_wallppfx"

DO_FOCUS="0"
FOCUS_DIR="Focussing"

DO_ALIGN="0"
ALIGN_DIR="Alignment"

# for i in nu nubar; do
for i in nu; do

  if [ "${DO_PPFX}" == "1" ]; then
    for dummy in dummy; do # So that continue does what you expect
    if [ -e /pnfs/dune/persistent/users/${USER}/${PPFX_OUTPUT_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i} ]; then
      echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/${USER}/${PPFX_OUTPUT_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i}"
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
      -p ${PPFX_OUTPUT_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i} \
      --number-of-jobs 3500 ${PPFX_ARG} \
      --expected-disk 1GB \
      --expected-mem 1.9GB \
      --expected-walltime 4h \
      --jobname BeamSim_${i}_PPFX
  done # end of dummy loop
  fi

  if [ "${DO_PPFX_COMPONENT_VARIATIONS}" == "1" ]; then
    for dummy in dummy; do # So that continue does what you expect

    if [ -e /pnfs/dune/persistent/users/${USER}/${PPFX_COMP_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i} ]; then
      echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/${USER}/${PPFX_COMP_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i}"
      continue
    fi

    if [ ! -e ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode.mac ]; then
      echo "[INFO]: No such g4lbnf macro OptimizedEngineeredNov2017Review_${i}mode.mac Not running."
      continue
    fi

    ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
      -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode.mac \
      -p ${PPFX_COMP_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i} \
      --number-of-jobs 2500 --PPFX --Use-All-PPFX-Weights \
      --expected-disk 2GB \
      --expected-mem 2.5GB \
      --expected-walltime 2h \
      --jobname BeamSim_${i}_PPFXAllW
  done #end of dummy loop
  fi

  #focussing
  if [ "${DO_FOCUS}" == "1" ]; then
    for j in p1 m1; do
      for k in WL HC DPR TargetDensity BeamSigma BeamOffsetX BeamTheta BeamThetaPhi; do

      if [ -e /pnfs/dune/persistent/users/${USER}/${FOCUS_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i} ]; then
        echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/${USER}/${FOCUS_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}"
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
       --expected-walltime 4h \
       --jobname BeamSim_${i}_${k}${j}
     done
   done
 fi

#alignment
  if [ "${DO_ALIGN}" == "1" ]; then
    for j in Horn1 Horn2; do
      for k in X Y XNeg X3mm XNeg3mm; do

      if [ -e /pnfs/dune/persistent/users/${USER}/${ALIGN_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i} ]; then
        echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/${USER}/${ALIGN_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i}"
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
       --expected-walltime 4h \
       --jobname BeamSim_${i}_${j}${k}Shift
      done
    done
  fi


#alignment
  if [ "${DO_HIGHERHC}" == "1" ]; then
    for j in 200 250; do

      if [ -e /pnfs/dune/persistent/users/${USER}/${HIGHERHC_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/HC_${j}/${i} ]; then
        echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/${USER}/${HIGHERHC_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/HC_${j}/${i}"
        continue
      fi

      if [ ! -e ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_HC${j}kA_${i}mode.mac ]; then
        echo "[INFO]: No such g4lbnf macro OptimizedEngineeredNov2017Review_HC${j}kA_${i}mode.mac Not running."
        continue
      fi

    EMEM=1GB
    PPFX_ARG=""
    if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
      PPFX_ARG="--PPFX"
      EMEM=1.9GB
    fi

      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
       -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_HC${j}kA_${i}mode.mac \
       -p ${HIGHERHC_DIR}/v3r5p4/QGSP_BERT/QGSP_BERT/OptimizedEngineeredNov2017Review/HC_${j}/${i} \
       --number-of-jobs 5000 \
       --expected-disk 1GB \
       --expected-mem ${EMEM} \
       --expected-walltime 4h \
       --jobname BeamSim_${i}_HC_${j} ${PPFX_ARG}
    done

    for j in 295.5 298 300.5 303 305.5 308 310.5 313 323 333 343; do

      if [ -e /pnfs/dune/persistent/users/${USER}/${HIGHERHC_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/HC_${j}/${i} ]; then
        echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/${USER}/${HIGHERHC_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/HC_${j}/${i}"
        continue
      fi

      if [ ! -e ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_HC${j}kA_${i}mode.mac ]; then
        echo "[INFO]: No such g4lbnf macro OptimizedEngineeredNov2017Review_HC${j}kA_${i}mode.mac Not running."
        continue
      fi

    EMEM=1GB
    PPFX_ARG=""
    if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
      PPFX_ARG="--PPFX"
      EMEM=1.9GB
    fi

      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
       -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_HC${j}kA_${i}mode.mac \
       -p ${HIGHERHC_DIR}/v3r5p4/QGSP_BERT/QGSP_BERT/OptimizedEngineeredNov2017Review/HC_${j}/${i} \
       --number-of-jobs 2000 \
       --expected-disk 1GB \
       --expected-mem ${EMEM} \
       --expected-walltime 4h \
       --jobname BeamSim_${i}_HC_${j} ${PPFX_ARG}
    done

  fi
done
