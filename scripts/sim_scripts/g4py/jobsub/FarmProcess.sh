#!/bin/bash

if [ -z ${DUNEPRISMTOOLSROOT} ]; then
  echo "[ERROR]: DUNEPRISM tools is not sourced, please set up the environment and try again."
  exit 1
fi

CONFIG=${DUNEPRISMTOOLSROOT}/configs/run_plans/RunPlan.39mLAr.4mx3mx5mActive.xml
CONDENSEDDIR=""
NMAXJOBS=""
OUTPUTDIR=""
ENVSETUPSCRIPT="${DUNEPRISMTOOLSROOT}/setup.sh"
FORCE="0"
POTPERFILE=""
SKIP="0"

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

      -R|--runplan-config)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      CONFIG="$2"
      echo "[OPT]: Using ${CONFIG} as the runplan configuration file."
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

      -S|--skip-first-N)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      SKIP="$2"
      echo "[OPT]: Skipping the first: \"${SKIP}\" files."
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

      -P|--POT-per-file)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      POTPERFILE="$2"
      echo "[OPT]: Using: \"${POTPERFILE}\" POT per file."
      shift # past argument
      ;;

      -f|--force)

      FORCE="1"
      echo "[OPT]: Will force overwrite of output."
      shift # past argument
      ;;

      -?|--help)

      echo "[RUNLIKE] ${SCRIPTNAME}"
      echo -e "\t-C|--condensed-input-dir   : Re-process condensed files found here."
      echo -e "\t-R|--runplan-config        : Use the specified runplan config xml."
      echo -e "\t-N|--NMAXJobs              : Maximum number of jobs to submit."
      echo -e "\t-S|--skip-first-N          : Number of files to skip when building joblist."
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
NINPUTSAFTERSKIP=$(( ${NINPUTS} - ${SKIP} ))

cat ${IPFL} | tail -${NINPUTSAFTERSKIP} > short.list
mv short.list ${IPFL}

NJOBSTORUN=$(( ${NINPUTS} - ${SKIP} ))
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

echo "qsub ${DUNEPRISMTOOLSROOT}/scripts/g4py/Process.sh -t 1-${NJOBSTORUN} -N DP_Proc \
  -v INPUT_FILE_LIST=${IPFL},RUNPLAN_CONFIG=${CONFIG},PROCESSED_OUTPUT_DIR=${PDIR},ENVSETUPSCRIPT=${ENVSETUPSCRIPT} \
  -o process.sub.log"

qsub ${DUNEPRISMTOOLSROOT}/scripts/g4py/Process.sh -t 1-${NJOBSTORUN} -N DP_Proc \
  -v INPUT_FILE_LIST=${IPFL},RUNPLAN_CONFIG=${CONFIG},PROCESSED_OUTPUT_DIR=${PDIR},ENVSETUPSCRIPT=${ENVSETUPSCRIPT},POTPERFILE=${POTPERFILE} \
  -o process.sub.log
