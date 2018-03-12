#!/bin/bash

if [ -z ${DUNEPRISMTOOLSROOT} ]; then
  echo "[ERROR]: DUNEPRISM tools is not sourced, please set up the environment and try again."
  exit 1
fi

CONFIG=${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml
CONDENSEDDIR=""
NMAXJOBS=""
OUTPUTDIR=""
ENVSETUPSCRIPT="${DUNEPRISMTOOLSROOT}/setup.sh"
FORCE="0"

while [[ ${#} -gt 0 ]]; do

  key="$1"
  case $key in

      -C|--condensed-input-dir)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      CONDENSEDDIR="$2"
      echo "[OPT]: Looking for condensed input files in ${CONDENSEDDIR}"
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
      shift # past argument
      ;;

      -?|--help)

      echo "[RUNLIKE] ${SCRIPTNAME}"
      echo -e "\t-G|--g4py-input-dir        : Build process list from all g4py files found here."
      echo -e "\t-R|--rootracker-input-dir  : Look here for associated rootracker files."
      echo -e "\t-N|--NMAXJobs              : Maximum number of jobs to submit."
      echo -e "\t-o|--output-dir            : Write output to <-o>/Condensed.<date> and <-o>/Processed.<date>."
      echo -e "\t-E|--env-script            : An environment setup script to source on the processing node {default: \${DUNEPRISMTOOLSROOT}/setup.sh}"
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

if [ -z ${CONDENSEDDIR} ] ; then
  echo "[ERROR]: Script must be passed both -C option. "
  exit 1
fi

if [ ! -e ${CONDENSEDDIR} ]; then
  echo "[ERROR]: Invalid -C directory \"${CONDENSEDDIR}\". "
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

source ${ENVSETUPSCRIPT}

if [ -z ${NMAXJOBS} ]; then
  ls ${CONDENSEDDIR}/*Condensed.root > condensed.process.list
else
  ls ${CONDENSEDDIR}/*Condensed.root | head -${NMAXJOBS} > condensed.process.list
fi

IPFL=$(readlink -f condensed.process.list)

NINPUTS=$(cat ${IPFL} | wc -l)

NJOBSTORUN=${NINPUTS}
if [ ! -z ${NMAXJOBS} ]; then
  NJOBSTORUN=$(python -c "print min(${NMAXJOBS},${NJOBSTORUN})")
fi

echo "[INFO]: Running ${NJOBSTORUN} job(s)."

PDIR="${OUTPUTDIR}/Processed.$(date +%Y-%m-%d)"

if [ -e ${PDIR} ]; then
  if [ ${FORCE} == 1 ]; then
    rm -r ${PDIR}
  else
    echo "[ERROR]: Output dir ${PDIR} already exists, will not overrwrite."
    exit 1
  fi
fi

mkdir -p ${PDIR}

if [ ! -e ${PDIR} ]; then
  echo "[ERROR]: Failed to make output dir: ${PDIR}"
  exit 1
fi

echo "qsub ${DUNEPRISMTOOLSROOT}/scripts/Process.sh -t 1-${NJOBSTORUN} -N DP_Proc \
  -v INPUT_FILE_LIST=${IPFL},RUNPLAN_CONFIG=${CONFIG},PROCESSED_OUTPUT_DIR=${PDIR},ENVSETUPSCRIPT=${ENVSETUPSCRIPT} \
  -o process.sub.log"

qsub ${DUNEPRISMTOOLSROOT}/scripts/Process.sh -t 1-${NJOBSTORUN} -N DP_Proc \
  -v INPUT_FILE_LIST=${IPFL},RUNPLAN_CONFIG=${CONFIG},PROCESSED_OUTPUT_DIR=${PDIR},ENVSETUPSCRIPT=${ENVSETUPSCRIPT} \
  -o process.sub.log
