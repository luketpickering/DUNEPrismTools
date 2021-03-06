# !/bin/bash

PRED_TYPE="uncertbin"
NDAXIS="offaxis"
declare -A FD_FHICL_ARR
declare -A ND_FHICL_ARR
FD_FHICL_ARR["nu"]="-FF flux/build_FD_fluxes_numode_${PRED_TYPE}_onaxis.fcl"
ND_FHICL_ARR["nu"]="-FN flux/build_ND_fluxes_numode_${PRED_TYPE}_${NDAXIS}.fcl"
FD_FHICL_ARR["nubar"]="-FF flux/build_FD_fluxes_nubarmode_${PRED_TYPE}_onaxis.fcl"
ND_FHICL_ARR["nubar"]="-FN flux/build_ND_fluxes_nubarmode_${PRED_TYPE}_${NDAXIS}.fcl"

PRED_DIR="uncertbin_offaxis"

FORCEOVERWRITE="true"

#coarsebin
# EDISK="2GB"
# EDISK_PPFX="2GB"
# EDISK_PPFX_COMP="2GB"

# # 3 files takes 30 s
# ETIME="30m"
# # 3 files takes 1 min
# ETIME_PPFX="2h"
# # 3 files takes 6 mins
# ETIME_PPFX_COMP="2h"

#finebin
EDISK="2GB"
EDISK_PPFX="2GB"
EDISK_PPFX_COMP="2GB"

# 3 files takes 1 min
ETIME="30m"
# 3 files takes 1 hour
ETIME_PPFX="2h"

EMEM="1GB"
EMEM_PPFX="1.9GB"
EMEM_PPFX_COMP="1.9GB"

declare -A PDG_ONLY
declare -A PDG_ONLY_PPFX
declare -A PDG_ONLY_PPFX_COMP

PDG_ONLY["nu"]="0"
PDG_ONLY_PPFX["nu"]="0"
PDG_ONLY_PPFX_COMP["nu"]="14"
PDG_ONLY["nubar"]="0"
PDG_ONLY_PPFX["nubar"]="0"
PDG_ONLY_PPFX_COMP["nubar"]="-14"

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

DO_HC280KA="1"
HC280KA_DIR="280KA"
if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
  HC280KA_DIR="280KA_wppfx"
fi
#coarsebin
NINPUTSPERJOB="50"
NMAXCONC="2"

NINPUTSPERJOB_PPFX="30"
NMAXCONC_PPFX="10"

NINPUTSPERJOB_PPFX_COMP="5"
NMAXCONC_PPFX_COMP="20"

#finebin
# NINPUTSPERJOB="20"
# NMAXCONC="5"

# NINPUTSPERJOB_PPFX="5"
# NMAXCONC_PPFX="20"

#finebin onaxis
# NINPUTSPERJOB="500"
# NMAXCONC="1"

# NINPUTSPERJOB_PPFX="250"
# NMAXCONC_PPFX="2"

NMAXJOBS="10000"
# NTOTALJOBS="1000000"

#BEAM_MODES="nu"
BEAM_MODES="nu nubar"

NJOBS="0"

for NUMODE in ${BEAM_MODES}; do

  for dummy in dummy; do

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
        --expected-walltime ${ETIME_PPFX} --expected-disk ${EDISK_PPFX} \
        --expected-mem ${EMEM_PPFX} ${FD_FHICL} ${ND_FHICL} \
        -p ${PPFX_DIR}/DUNEPrismFluxes/__DET___${NUMODE}/${PRED_DIR} \
        -i /pnfs/dune/persistent/users/${USER}/${PPFX_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${NUMODE}/dk2nulite \
        -n ${NINPUTSPERJOB_PPFX} --maxConcurrent ${NMAXCONC_PPFX}\
        --N-max-jobs ${NMAXJOBS} ${PPFX_ARG}\
        --only-pdg ${PDG_ONLY_PPFX[${NUMODE}]} -f
      NJOBS=$(( NJOBS + NMAXJOBS ))
      echo "[INFO]: Done ${NJOBS}/${NTOTALJOBS}."
      if [ ${NJOBS} -ge ${NTOTALJOBS} ]; then
        echo "[INFO]: Stopping after ${NJOBS} JOBS"
        exit
      fi
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
        -n ${NINPUTSPERJOB_PPFX_COMP} --maxConcurrent ${NMAXCONC_PPFX_COMP}\
          --N-max-jobs ${NMAXJOBS} ${PPFX_ARG}\
        --only-pdg ${PDG_ONLY_PPFX_COMP[${NUMODE}]} -f

      NJOBS=$(( NJOBS + NMAXJOBS ))
      echo "[INFO]: Done ${NJOBS}."
      if [ ${NJOBS} -ge ${NTOTALJOBS} ]; then
        echo "[INFO]: Stopping after ${NJOBS} JOBS"
        exit
      fi
    fi
  done

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
           -n ${NINPUTSPERJOB} --maxConcurrent ${NMAXCONC}\
           --N-max-jobs ${NMAXJOBS} \
           --only-pdg ${PDG_ONLY[${NUMODE}]} -f

        NJOBS=$(( NJOBS + NMAXJOBS ))
        echo "[INFO]: Done ${NJOBS}."
        if [ ${NJOBS} -ge ${NTOTALJOBS} ]; then
          echo "[INFO]: Stopping after ${NJOBS} JOBS"
          exit
        fi
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
           -n ${NINPUTSPERJOB} --maxConcurrent ${NMAXCONC}\
            --N-max-jobs ${NMAXJOBS} \
           --only-pdg ${PDG_ONLY[${NUMODE}]} -f

        NJOBS=$(( NJOBS + NMAXJOBS ))
        echo "[INFO]: Done ${NJOBS}."
        if [ ${NJOBS} -ge ${NTOTALJOBS} ]; then
          echo "[INFO]: Stopping after ${NJOBS} JOBS"
          exit
        fi
      done
    done
  fi

  #Alignment
  if [ "${DO_HC280KA}" == "1" ]; then
    if [ ! -e /pnfs/dune/persistent/users/${USER}/${HC280KA_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${NUMODE}/dk2nulite ]; then
      echo "[INFO]: No input directory, skipping."
      continue
    fi

    if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/${USER}/${HC280KA_DIR}/DUNEPrismFluxes/ND_${NUMODE}/${PRED_DIR} ]; then
      echo "[INFO]: Already have ND_${NUMODE}/280kA not reprocessing."
      continue
    fi


    ETIME_280KA=${ETIME}
    EDISK_280KA=${EDISK}
    EMEM_280KA=${EMEM}
    NMAXCONC_280KA=${NMAXCONC}
    if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
      ETIME_280KA=${ETIME_PPFX}
      EDISK_280KA=${EDISK_PPFX}
      EMEM_280KA=${EMEM_PPFX}
      NMAXCONC_280KA=${NMAXCONC_PPFX}
    fi

    PPFX_ARG=""
    if [ "${DO_PPFX_VARIATIONS}" == "1" ]; then
      PPFX_ARG="--NPPFXU ${NPPFXUNIVERSES}"
    fi

    ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
       --expected-walltime ${ETIME_280KA} --expected-disk ${EDISK_280KA} \
       --expected-mem ${EMEM_280KA} ${FD_FHICL} ${ND_FHICL} \
       -p ${HC280KA_DIR}/DUNEPrismFluxes/__DET___${NUMODE}/${PRED_DIR} \
       -i /pnfs/dune/persistent/users/${USER}/${HC280KA_DIR}/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${NUMODE}/dk2nulite \
       -n ${NINPUTSPERJOB} --maxConcurrent ${NMAXCONC_280KA}\
        --N-max-jobs ${NMAXJOBS} ${PPFX_ARG}\
       --only-pdg ${PDG_ONLY[${NUMODE}]} -f

    NJOBS=$(( NJOBS + NMAXJOBS ))
    echo "[INFO]: Done ${NJOBS}."
    if [ ${NJOBS} -ge ${NTOTALJOBS} ]; then
      echo "[INFO]: Stopping after ${NJOBS} JOBS"
      exit
    fi
  fi

done
