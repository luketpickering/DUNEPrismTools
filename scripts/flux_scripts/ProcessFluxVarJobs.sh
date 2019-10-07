# !/bin/bash

# PRED_DIR="finebin"
PRED_DIR="uncertbin"
declare -A FD_FHICL_ARR
declare -A ND_FHICL_ARR
FD_FHICL_ARR["nu"]="-FF flux/build_FD_fluxes_numode_${PRED_DIR}_onaxis.fcl"
ND_FHICL_ARR["nu"]="-FN flux/build_ND_fluxes_numode_${PRED_DIR}_offaxis.fcl"
FD_FHICL_ARR["nubar"]="-FF flux/build_FD_fluxes_nubarmode_${PRED_DIR}_onaxis.fcl"
ND_FHICL_ARR["nubar"]="-FN flux/build_ND_fluxes_nubarmode_${PRED_DIR}_offaxis.fcl"

FORCEOVERWRITE="true"

#finebin
# EDISK="4GB"
# EDISK_PPFX_COMP="6GB"
# ETIME_PPFX="4h"
# ETIME_PPFX_COMP="10h"
# ETIME="2h"
# EMEM="512MB"
# EMEM_PPFX="1GB"
# EMEM_PPFX_COMP="2GB"

#coarsebin
EDISK="4GB"
EDISK_PPFX_COMP="8GB"
ETIME_PPFX="2h"
ETIME_PPFX_COMP="12h"
ETIME="2h"
EMEM="512MB"
EMEM_PPFX="1GB"
EMEM_PPFX_COMP="2GB"

DO_PPFX="0"
NPPFXUNIVERSES="100"
PPFX_DIR="nominal_5E8POT_wppfx"
DO_PPFX_VARIATIONS="1"

DO_PPFX_COMPONENT_VARIATIONS="1"
PPFX_COMP_DIR="nominal_2.5E8POT_wallppfx"

DO_HIGHERHC="1"
HIGHERHC_DIR="HigherHC_2.5E8POT"
if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
  HIGHERHC_DIR="HigherHC_2.5E8POT_wppfx"
fi

DO_FOCUS="0"
FOCUS_DIR="Focussing"

DO_ALIGN="0"
ALIGN_DIR="Alignment"

NINPUTSPERJOB="1"
NMAXJOBS="1"

for NUMODE in nu nubar; do

  FD_FHICL=${FD_FHICL_ARR[${NUMODE}]}
  ND_FHICL=${ND_FHICL_ARR[${NUMODE}]}

  if [ "${DO_PPFX}" == "1" ]; then
    if [ ! -e /pnfs/dune/persistent/users/${USER}/${PPFX_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${NUMODE}/dk2nulite ]; then
      echo "[INFO]: No input directory for ${NUMODE} wppfx, skipping."
      continue
    fi

    if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/${USER}/${PPFX_DIR}/DUNEPrismFluxes/ND_${NUMODE}/${PRED_DIR} ]; then
      echo "[INFO]: Already have ${NUMODE} wppfx not reprocessing."
      continue
    fi

    PPFX_ARG=""
    if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
      PPFX_ARG="--NPPFXU ${NPPFXUNIVERSES}"
    fi

  ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
     --expected-walltime ${ETIME_PPFX} --expected-disk ${EDISK} \
     --expected-mem ${EMEM_PPFX} ${FD_FHICL} ${ND_FHICL} \
     -p ${PPFX_DIR}/DUNEPrismFluxes/__DET___${NUMODE}/${PRED_DIR} \
     -i /pnfs/dune/persistent/users/${USER}/${PPFX_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${NUMODE}/dk2nulite \
     -n ${NINPUTSPERJOB} --N-max-jobs ${NMAXJOBS} ${PPFX_ARG} -f
  fi

  if [ "${DO_PPFX_COMPONENT_VARIATIONS}" == "1" ]; then
    if [ ! -e /pnfs/dune/persistent/users/${USER}/${PPFX_COMP_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${NUMODE}/dk2nulite ]; then
      echo "[INFO]: No input directory for ${NUMODE} wppfx, skipping."
      continue
    fi

    if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/${USER}/${PPFX_COMP_DIR}/DUNEPrismFluxes/ND_${NUMODE}/${PRED_DIR} ]; then
      echo "[INFO]: Already have ${NUMODE} wppfx not reprocessing."
      continue
    fi

    PPFX_ARG=""
    if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
      PPFX_ARG="--NPPFXU ${NPPFXUNIVERSES} --PPFX-Components"
    fi

  ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
     --expected-walltime ${ETIME_PPFX_COMP} --expected-disk ${EDISK_PPFX_COMP} \
     --expected-mem ${EMEM_PPFX_COMP} ${FD_FHICL} ${ND_FHICL} \
     -p ${PPFX_COMP_DIR}/DUNEPrismFluxes/__DET___${NUMODE}/${PRED_DIR} \
     -i /pnfs/dune/persistent/users/${USER}/${PPFX_COMP_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${NUMODE}/dk2nulite \
     -n ${NINPUTSPERJOB} --N-max-jobs ${NMAXJOBS} ${PPFX_ARG} -f
  fi

  #With focussing
  if [ "${DO_FOCUS}" == "1" ]; then
    for SIGSHIFT in p1 m1; do
      for VARIATION in WL HC DPR TargetDensity BeamSigma BeamOffsetX BeamTheta BeamThetaPhi; do

        if [ ! -e /pnfs/dune/persistent/users/${USER}/${FOCUS_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${VARIATION}${SIGSHIFT}/${NUMODE}/dk2nulite ]; then
          echo "[INFO]: No input directory, skipping."
          continue
        fi

        if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/${USER}/${FOCUS_DIR}/DUNEPrismFluxes/ND_${NUMODE}/${VARIATION}${SIGSHIFT}/${PRED_DIR} ]; then
          echo "[INFO]: Already have ND_${NUMODE}/${VARIATION}${SIGSHIFT} not reprocessing."
          continue
        fi

        ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
           --expected-walltime ${ETIME} --expected-disk ${EDISK} \
           --expected-mem ${EMEM} ${FD_FHICL} ${ND_FHICL} \
           -p ${FOCUS_DIR}/DUNEPrismFluxes/__DET___${NUMODE}/${VARIATION}${SIGSHIFT}/${PRED_DIR} \
           -i /pnfs/dune/persistent/users/${USER}/${FOCUS_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${VARIATION}${SIGSHIFT}/${NUMODE}/dk2nulite \
           -n ${NINPUTSPERJOB} --N-max-jobs ${NMAXJOBS} -f
      done
    done
  fi

  #Alignment
  if [ "${DO_ALIGN}" == "1" ]; then
    for HORN in Horn1 Horn2; do
      for VARIATION in X Y XNeg X3mm XNeg3mm; do

        if [ ! -e /pnfs/dune/persistent/users/${USER}/${ALIGN_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${HORN}${VARIATION}/${NUMODE}/dk2nulite ]; then
          echo "[INFO]: No input directory, skipping."
          continue
        fi

        if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/${USER}/${ALIGN_DIR}/DUNEPrismFluxes/ND_${NUMODE}/${HORN}${VARIATION}/${PRED_DIR} ]; then
          echo "[INFO]: Already have ND_${NUMODE}/${VARIATION}${HORN} not reprocessing."
          continue
        fi

        ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
           --expected-walltime ${ETIME} --expected-disk ${EDISK} \
           --expected-mem ${EMEM} ${FD_FHICL} ${ND_FHICL} \
           -p ${ALIGN_DIR}/DUNEPrismFluxes/__DET___${NUMODE}/${HORN}${VARIATION}/${PRED_DIR} \
           -i /pnfs/dune/persistent/users/${USER}/${ALIGN_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${HORN}${VARIATION}/${NUMODE}/dk2nulite \
           -n ${NINPUTSPERJOB} --N-max-jobs ${NMAXJOBS} -f
      done
    done
  fi


    #Alignment
  if [ "${DO_HIGHERHC}" == "1" ]; then
    for CURR in 303 313 323 333 343; do
        if [ ! -e /pnfs/dune/persistent/users/${USER}/${HIGHERHC_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/HC_${CURR}/${NUMODE}/dk2nulite ]; then
          echo "[INFO]: No input directory, skipping."
          continue
        fi

        if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/${USER}/${HIGHERHC_DIR}/DUNEPrismFluxes/ND_${NUMODE}/HC_${CURR}/${PRED_DIR} ]; then
          echo "[INFO]: Already have ND_${NUMODE}/${VARIATION}${HORN} not reprocessing."
          continue
        fi

        HHC_ETIME=${ETIME}
        HHC_EDISK=${EDISK}
        HHC_EMEM=${EMEM}

        PPFX_ARG=""
        if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
          PPFX_ARG="--NPPFXU ${NPPFXUNIVERSES}"
          HHC_ETIME=${ETIME_PPFX}
          HHC_EDISK=${EDISK_PPFX}
          HHC_EMEM=${EMEM_PPFX}
        fi

        ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
           --expected-walltime ${HHC_ETIME} --expected-disk ${HHC_EDISK} \
           --expected-mem ${HHC_EMEM} ${FD_FHICL} ${ND_FHICL} \
           -p ${HIGHERHC_DIR}/DUNEPrismFluxes/__DET___${NUMODE}/HC_${CURR}/${PRED_DIR} \
           -i /pnfs/dune/persistent/users/${USER}/${ALIGN_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/HC_${CURR}/${NUMODE}/dk2nulite \
           -n ${NINPUTSPERJOB} --N-max-jobs ${NMAXJOBS} -f ${PPFX_ARG}
      done
    done
  fi

done
