#!/bin/bash

G4PYSCRATCHDIR=""
ROOTRACKERDIR=""
OUPFILE="process.list"
FORCE="0"
NMAX=""
SKIP="0"

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

      -f|--force)


      FORCE="1"
      echo "[OPT]: Skipping g4py files without corresponding rootracker files."
      ;;

      -o|--output-file)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      OUPFILE="$2"
      echo "[OPT]: Writing output list to ${OUPFILE}"
      shift # past argument
      ;;

      -N|--nmax-files)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      NMAX="$2"
      echo "[OPT]: Reading a maximum of ${NMAX} file"
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

      -?|--help)

      echo "[RUNLIKE] ${SCRIPTNAME}"
      echo -e "\t-G|--g4py-input-dir        : Build process list from all g4py files found here."
      echo -e "\t-R|--rootracker-input-dir  : Look here for associated rootracker files."
      echo -e "\t-f|--force                 : Will not stop for missing rootracker files."
      echo -e "\t-o|--output-file           : output file name {default \"process.list\"}."
      echo -e "\t-N|--nmax-files            : Maximum number of files to find."
      echo -e "\t-S|--skip-first-N          : Number of files to skip when building joblist."
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

if [ -e ${OUPFILE} ]; then
  rm ${OUPFILE}
fi
touch ${OUPFILE}

FILES=$(ls ${G4PYSCRATCHDIR}/*.root)
NFOUND=0

if [ ! -z ${NMAX} ]; then
  NMAX=$(( ${NMAX} + ${SKIP} ))
fi

for F in ${FILES}; do
  echo "[INFO]:===="
  echo "[INFO]: File = ${F}"
  echo "[INFO]: Filename = ${F##*/}"
  PROCID=$( echo ${F##*/} | awk '{split($0,a,"."); print a[1]}' )
  CLUSTID=$( echo ${F##*/} | awk '{split($0,a,"."); print a[2]}' )

  if [[ ${NFOUND} -ge ${SKIP} ]]; then

    echo "[INFO]: Looking for rootracker file g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredSept2017Review_*.rootracker.${PROCID}.${CLUSTID}.root"

    if ! ls ${ROOTRACKERDIR}/g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredSept2017Review_*.rootracker.${PROCID}.${CLUSTID}.root &> /dev/null; then
      if [ ${FORCE} == "0" ]; then
        echo "[ERROR]: Failed to find matching rootracker file, cannot continue. Pass -f to skip input instead."
        exit 1
      else
        continue
      fi
    fi
    RTFILE=$(ls ${ROOTRACKERDIR}/g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredSept2017Review_*.rootracker.${PROCID}.${CLUSTID}.root)

    echo "[INFO]:===="

    echo "${F} ${RTFILE}" >> ${OUPFILE}
  else
    echo "[INFO]: Skipping file [${NFOUND}] ${F##*/}"
  fi

  (( NFOUND++ ))
  if [ ! -z ${NMAX} ]; then
    if [ ${NFOUND} == ${NMAX} ]; then
      echo "[INFO]: Found ${NFOUND} files, stopping search."
      break
    fi
  fi
done
