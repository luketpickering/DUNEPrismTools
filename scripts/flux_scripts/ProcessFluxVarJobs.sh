# !/bin/bash

PRED_DIR="uncertbin_onaxis"
declare -A FD_FHICL_ARR
declare -A ND_FHICL_ARR
FD_FHICL_ARR["nu"]="-FF flux/build_FD_fluxes_numode_uncertbin_onaxis.fcl"
ND_FHICL_ARR["nu"]="-FN flux/build_ND_fluxes_numode_uncertbin_onaxis.fcl"
FD_FHICL_ARR["nubar"]="-FF flux/build_FD_fluxes_nubarmode_uncertbin_onaxis.fcl"
ND_FHICL_ARR["nubar"]="-FN flux/build_ND_fluxes_nubarmode_uncertbin_onaxis.fcl"

FORCEOVERWRITE="false"

EDISK="4GB"
ETIME_PPFX="6h"
ETIME="2h"

DO_PPFX="1"
NPPFXUNIVERSES="100"
PPFX_DIR="nominal_5E8POT_wppfx"
DO_PPFX_VARIATIONS="1"

DO_PPFX_COMPONENT_VARIATIONS="0"
PPFX_COMP_DIR="nominal_2.5E8POT_wallppfx"

DO_FOCUS="1"
FOCUS_DIR="Focussing"

DO_ALIGN="1"
ALIGN_DIR="Alignment"

NINPUTSPERJOB="15"

for NUMODE in nu nubar; do

  FD_FHICL=${FD_FHICL_ARR[${NUMODE}]}
  ND_FHICL=${ND_FHICL_ARR[${NUMODE}]}

  if [ "${DO_PPFX}" == "1" ]; then
    if [ ! -e /pnfs/dune/persistent/users/picker24/${PPFX_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${NUMODE}/dk2nulite ]; then
      echo "[INFO]: No input directory for ${NUMODE} wppfx, skipping."
      continue
    fi

    if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/picker24/${PPFX_DIR}/DUNEPrismFluxes/ND_${NUMODE}/${PRED_DIR} ]; then
      echo "[INFO]: Already have ${NUMODE} wppfx not reprocessing."
      continue
    fi

    PPFX_ARG=""
    if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
      PPFX_ARG="--NPPFXU ${NPPFXUNIVERSES}"
    fi

  ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
     --expected-walltime ${ETIME_PPFX} --expected-disk ${EDISK} \
     --expected-mem 512MB ${FD_FHICL} ${ND_FHICL} \
     -p ${PPFX_DIR}/DUNEPrismFluxes/__DET___${NUMODE}/${PRED_DIR} \
     -i /pnfs/dune/persistent/users/picker24/${PPFX_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${NUMODE}/dk2nulite \
     -n ${NINPUTSPERJOB} ${PPFX_ARG} -f
  fi

  if [ "${DO_PPFX_COMPONENT_VARIATIONS}" == "1" ]; then
    if [ ! -e /pnfs/dune/persistent/users/picker24/${PPFX_COMP_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${NUMODE}/dk2nulite ]; then
      echo "[INFO]: No input directory for ${NUMODE} wppfx, skipping."
      continue
    fi

    if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/picker24/${PPFX_COMP_DIR}/DUNEPrismFluxes/ND_${NUMODE}/${PRED_DIR} ]; then
      echo "[INFO]: Already have ${NUMODE} wppfx not reprocessing."
      continue
    fi

    PPFX_ARG=""
    if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
      PPFX_ARG="--NPPFXU ${NPPFXUNIVERSES}"
    fi

  ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
     --expected-walltime ${ETIME_PPFX} --expected-disk ${EDISK} \
     --expected-mem 512MB ${FD_FHICL} ${ND_FHICL} \
     -p ${PPFX_COMP_DIR}/DUNEPrismFluxes/__DET___${NUMODE}/${PRED_DIR} \
     -i /pnfs/dune/persistent/users/picker24/${PPFX_COMP_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${NUMODE}/dk2nulite \
     -n ${NINPUTSPERJOB} ${PPFX_ARG} -f
  fi

  #With focussing
  if [ "${DO_FOCUS}" == "1" ]; then
    for SIGSHIFT in p1 m1; do
      for VARIATION in WL HC DPR TargetDensity BeamSigma BeamOffsetX BeamTheta BeamThetaPhi; do

        if [ ! -e /pnfs/dune/persistent/users/picker24/${FOCUS_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${VARIATION}${SIGSHIFT}/${NUMODE}/dk2nulite ]; then
          echo "[INFO]: No input directory, skipping."
          continue
        fi

        if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/picker24/${FOCUS_DIR}/DUNEPrismFluxes/ND_${NUMODE}/${VARIATION}${SIGSHIFT}/${PRED_DIR} ]; then
          echo "[INFO]: Already have ND_${NUMODE}/${VARIATION}${SIGSHIFT} not reprocessing."
          continue
        fi

        ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
           --expected-walltime ${ETIME} --expected-disk ${EDISK} \
           --expected-mem 512MB ${FD_FHICL} ${ND_FHICL} \
           -p ${FOCUS_DIR}/DUNEPrismFluxes/__DET___${NUMODE}/${VARIATION}${SIGSHIFT}/${PRED_DIR} \
           -i /pnfs/dune/persistent/users/picker24/${FOCUS_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${VARIATION}${SIGSHIFT}/${NUMODE}/dk2nulite \
           -n ${NINPUTSPERJOB} -f
      done
    done
  fi

  #Alignment
  if [ "${DO_ALIGN}" == "1" ]; then
    for HORN in Horn1 Horn2; do
      for VARIATION in X Y XNeg X3mm XNeg3mm; do

        if [ ! -e /pnfs/dune/persistent/users/picker24/${ALIGN_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${HORN}${VARIATION}/${NUMODE}/dk2nulite ]; then
          echo "[INFO]: No input directory, skipping."
          continue
        fi

        if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/picker24/${ALIGN_DIR}/DUNEPrismFluxes/ND_${NUMODE}/${HORN}${VARIATION}/${PRED_DIR} ]; then
          echo "[INFO]: Already have ND_${NUMODE}/${VARIATION}${HORN} not reprocessing."
          continue
        fi

        ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
           --expected-walltime ${ETIME} --expected-disk ${EDISK} \
           --expected-mem 512MB ${FD_FHICL} ${ND_FHICL} \
           -p ${ALIGN_DIR}/DUNEPrismFluxes/__DET___${NUMODE}/${HORN}${VARIATION}/${PRED_DIR} \
           -i /pnfs/dune/persistent/users/picker24/${ALIGN_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${HORN}${VARIATION}/${NUMODE}/dk2nulite \
           -n ${NINPUTSPERJOB} -f
      done
    done
  fi
done
