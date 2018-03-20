#!/bin/bash

if [ -z ${DUNEPRISMTOOLSROOT} ]; then
  echo "[ERROR]: DUNEPRISM tools is not sourced, please set up the environment and try again."
  exit 1
fi

CONFIG=${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml
G4PYSCRATCHDIR=""
ROOTRACKERDIR=""
NMAXJOBS=""
OUTPUTDIR=""
ENVSETUPSCRIPT="${DUNEPRISMTOOLSROOT}/setup.sh"
FORCE="0"
CONDENSERCONFIG="-nx 400 -dmn -3800,-150,-250 -dmx 200,150,250 -fv 50,50,50 -nt 1000 -T 250"
POTPERFILE=""

while [[ ${#} -gt 0 ]]; do

  key="$1"
  case $key in

      -G|--g4py-input-dir)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      G4PYSCRATCHDIR="$2"
      echo "[OPT]: Reading g4py files from ${G4PYSCRATCHDIR}"
      shift # past argument
      ;;

      -R|--rootracker-input-dir)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      ROOTRACKERDIR="$2"
      echo "[OPT]: Looking for rootracker files from ${ROOTRACKERDIR}"
      shift # past argument
      ;;

      -C|--runplan-config)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      CONFIG="$2"
      echo "[OPT]: Will use ${CONFIG} as runplan configuration"
      shift # past argument
      ;;

      -F|--condenser-config)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      CONDENSERCONFIG="$2"
      echo "[OPT]: Will use ${CONDENSERCONFIG} as condenser configuration"
      shift # past argument
      ;;

      -N|--N-max-jobs)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      NMAXJOBS="$2"
      echo "[OPT]: Running a maximum of: \"${NMAXJOBS}\"."
      shift # past argument
      ;;

      -o|--output-dir)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      OUTPUTDIR="$2"
      echo "[OPT]: Writing output to: \"${OUTPUTDIR}\"."
      shift # past argument
      ;;

      -P|--POT-per-file)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      POTPERFILE="$2"
      echo "[OPT]: Using: \"${POTPERFILE}\" POT per file."
      shift # past argument
      ;;

      -E|--env-script)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      ENVSETUPSCRIPT="$2"
      echo "[OPT]: Using: \"${ENVSETUPSCRIPT}\" to set up environment."
      shift # past argument
      ;;

      -f|--force)

      FORCE="1"
      echo "[OPT]: Will force overwrite of output."
      ;;

      -?|--help)

      echo "[RUNLIKE] ${SCRIPTNAME}"
      echo -e "\t-G|--g4py-input-dir        : Build process list from all g4py files found here."
      echo -e "\t-R|--rootracker-input-dir  : Look here for associated rootracker files."
      echo -e "\t-C|--runplan-config        : Override default runplan configuration."
      echo -e "\t-N|--NMAXJobs              : Maximum number of jobs to submit."
      echo -e "\t-o|--output-dir            : Write output to <-o>/Condensed.<date> and <-o>/Processed.<date>."
      echo -e "\t-E|--env-script            : An environment setup script to source on the processing node {default: \${DUNEPRISMTOOLSROOT}/setup.sh}"
      echo -e "\t-P|--POT-per-file          : POT per processed file."
      echo -e "\t-f|--force                 : Will remove previous output directories if they exist"
      echo -e "\t-?|--help                  : Print this message."
      exit 0
      ;;

      *)
              # unknown option
      echo "Unknown option $1"
      exit 1
      ;;
  esac
  shift # past argument or value
done

if [ -z ${G4PYSCRATCHDIR} ] || [ -z ${ROOTRACKERDIR} ]; then
  echo "[ERROR]: Script must be passed both -G and -R options. "
  exit 1
fi

if [ ! -e ${G4PYSCRATCHDIR} ]; then
  echo "[ERROR]: Invalid -G directory \"${G4PYSCRATCHDIR}\". "
  exit 1
fi

if [ ! -e ${ROOTRACKERDIR} ]; then
  echo "[ERROR]: Invalid -R directory \"${ROOTRACKERDIR}\". "
  exit 1
fi

if [ -z ${OUTPUTDIR} ]; then
  echo "[ERROR]: No output directory passed, please supply a -o value."
  exit 1
fi

if [ -z ${ENVSETUPSCRIPT} ]; then
  echo "[ERROR]: No environment set up script passed."
  exit 1
fi

if [ ! -e ${CONFIG} ]; then
  echo "[ERROR]: Config script \"${CONFIG}\" doesn't exist."
  exit 1
fi

source ${ENVSETUPSCRIPT}

if [ -z ${NMAXJOBS} ]; then
  ${DUNEPRISMTOOLSROOT}/scripts/BuildInputsList_ProcessFromG4Ar.sh -G ${G4PYSCRATCHDIR} -R ${ROOTRACKERDIR} -o process.list
else
  ${DUNEPRISMTOOLSROOT}/scripts/BuildInputsList_ProcessFromG4Ar.sh -G ${G4PYSCRATCHDIR} -R ${ROOTRACKERDIR} -o process.list -N ${NMAXJOBS}
fi

IPFL=$(readlink -f process.list)

NINPUTS=$(cat ${IPFL} | wc -l)

NJOBSTORUN=${NINPUTS}
if [ ! -z ${NMAXJOBS} ]; then
  NJOBSTORUN=$(python -c "print min(${NMAXJOBS},${NJOBSTORUN})")
fi

CDIR="${OUTPUTDIR}/Condensed.$(date +%Y-%m-%d)"
PDIR="${OUTPUTDIR}/Processed.$(date +%Y-%m-%d)"

if [ -e ${CDIR} ]; then
  if [ ${FORCE} == 1 ]; then
    rm -r ${CDIR}
  else
    echo "[ERROR]: Output dir ${CDIR} already exists, will not overrwrite."
    exit 1
  fi
fi

if [ -e ${PDIR} ]; then
  if [ ${FORCE} == 1 ]; then
    rm -r ${PDIR}
  else
    echo "[ERROR]: Output dir ${PDIR} already exists, will not overrwrite."
    exit 1
  fi
fi

mkdir -p ${CDIR}
mkdir -p ${PDIR}

if [ ! -e ${CDIR} ]; then
  echo "[ERROR]: Failed to make output dir: ${CDIR}"
  exit 1
fi

if [ ! -e ${PDIR} ]; then
  echo "[ERROR]: Failed to make output dir: ${PDIR}"
  exit 1
fi

ENCCC=$(echo ${CONDENSERCONFIG} | tr " " "_" | tr "," "|")

echo "[INFO]: Encoded condenser config = \"${ENCCC}\""

qsub ${DUNEPRISMTOOLSROOT}/scripts/CondenseAndProcess.sh -t 1-${NJOBSTORUN} -N DP_CondAndProc \
  -v INPUT_FILE_LIST=${IPFL},RUNPLAN_CONFIG=${CONFIG},CONDENSED_OUTPUT_DIR=${CDIR},PROCESSED_OUTPUT_DIR=${PDIR},ENVSETUPSCRIPT=${ENVSETUPSCRIPT},CONDENSERCONFIG=${ENCCC},POTPERFILE=${POTPERFILE} \
  -o sub.log
