#!/bin/bash

if [ -z ${INPUT_TAR_FILE} ]; then
  echo "[ERROR]: Expected to recieve an input file."
  exit 1
fi

if [ ! -e dp_BuildFluxes ]; then
  echo "[ERROR]: Couldn't find dp_BuildFluxes"
  ls
  exit 2
fi

if [ ! -e inputs.list ]; then
  echo "[ERROR]: Couldn't find inputs.list"
  ls
  exit 3
fi

if [ -z $1 ]; then
  echo "[ERROR]: Expected to recieve a detector distance in cm as a script argument."
  exit 5
fi
DET_DIST_CM=${1}

if [ -z ${2} ]; then
  echo "[ERROR]: Couldn't find binning descriptor."
  exit 6
fi
BINNING_DESCRIPTOR=${2}

if [ -z ${3} ]; then
  echo "[ERROR]: Couldn't find Reuse parent descriptor."
  exit 6
fi
REUSEPARENTS=${3}

RUPARG=""
if [ "${REUSEPARENTS}" == "0" ]; then
  RUPARG="-P"
  echo "[INFO]: Not re-using decay parents."
else
  echo "[INFO]: Re-using decay parents."
fi

if [ -z ${4} ]; then
  echo "[ERROR]: Couldn't find species descriptor."
  exit 6
fi
ONLYSPECIES=${4}

SPECARG=""
if [ "${ONLYSPECIES}" != "0" ]; then
  SPECARG="-S ${ONLYSPECIES}"
  echo "[INFO]: Only building for species ${ONLYSPECIES}."
fi

if [ -z ${5} ]; then
  echo "[ERROR]: Couldn't find species descriptor."
  exit 6
fi
DK2NULITE=${5}

DK2NU_LITE_ARG=""
if [ "${DK2NULITE}" != "0" ]; then
  DK2NU_LITE_ARG="-L"
  echo "[INFO]: Expecting dk2nu_lite inputs."
fi

PNFS_PATH_APPEND=""
if [ ! -z ${6} ]; then
  echo "[INFO]: Appending path to pnfs outdir: ${3}"
  PNFS_PATH_APPEND=${6}
fi

NFILES_TO_READ=1
if [ ! -z ${7} ]; then
  NFILES_TO_READ=${7}
fi

if [ -z ${8} ]; then
  echo "[ERROR]: Expected to recieve a flux window descriptor."
  exit 1
fi
FLUX_WINDOW_DESCRIPTOR_ENC=${8}
FLUX_WINDOW_DESCRIPTOR=$(python -c "import urllib; print urllib.unquote('''${FLUX_WINDOW_DESCRIPTOR_ENC}''')")

echo "Unencoded flux window descriptor: ${FLUX_WINDOW_DESCRIPTOR}"

PPFX_ARG="${9}"
echo "[INFO]: PPFX Argument: \"${PPFX_ARG}\""

echo "[INFO]: Reading ${NFILES_TO_READ} files."

echo "[INFO]: JobID ${CLUSTER}, ArrayID ${PROCESS}"

(( ZERO_INDEX_LN = ${PROCESS} * ${NFILES_TO_READ} ))

(( LN = ${ZERO_INDEX_LN} + 1 ))
(( LAST_READ_LN = ${LN} + (${NFILES_TO_READ} - 1) ))

NTOT_FILES=$(cat inputs.list | wc -l)

echo "[INFO]: ZERO_INDEX_LN: ${ZERO_INDEX_LN}, LN: ${LN}, LAST_READ_LN: ${LAST_READ_LN}, NTOT_FILES: ${NTOT_FILES}"

if [ ${END_LN} > ${NTOT_FILES} ]; then
  ONF=${NFILES_TO_READ}
  (( NFILES_TO_READ = ${NFILES_TO_READ} - (${LAST_READ_LN} - ${NTOT_FILES}) ))
  LAST_READ_LN=${NTOT_FILES}
  echo "[INFO]: Updating files to read so that we don't fall off the end or process the same file twice: ${ONF} => ${NFILES_TO_READ}."
fi

INPUT_FILE_PATHS=( $( cat inputs.list | head -${LAST_READ_LN} | tail -${NFILES_TO_READ} ) )
NINPUT_FILES_FOUND=${#INPUT_FILE_PATHS[@]}

echo "[INFO]: cat inputs.list | head -${LN} | tail -${NFILES_TO_READ}"
echo "[INFO]: INPUT_FILE_PATHS: ${INPUT_FILE_PATHS[@]}"
echo "[INFO]: NINPUT_FILES_FOUND: ${NINPUT_FILES_FOUND}"

if [ -z ${NINPUT_FILES_FOUND} ]; then
  echo "[ERROR]: Found no input files with: cat inputs.list | head -${LN} | tail -${NFILES_TO_READ}"
  exit 6
fi

echo "Detector distance: ${DET_DIST_CM} "

printenv

set -x #start bash debugging at this point
echo "Start $(date)"
echo "Site:${GLIDEIN_ResourceName}"
echo "The worker node is " `hostname` "OS: " `uname -a`
echo "The user id is $(whoami)"
echo "The output of id is: $(id)"
set +x #stop bash debugging at this point

if [ -z ${GRID_USER} ]; then
  GRID_USER=$(basename $X509_USER_PROXY | cut -d "_" -f 2)
fi

if [ -z ${GRID_USER} ]; then
  echo "Failed to get GRID_USER."
  exit 9
fi

cp dp_BuildFluxes $_CONDOR_SCRATCH_DIR/

cd $_CONDOR_SCRATCH_DIR

echo "pwd is $(pwd)"
echo "------ls-------"
ls
echo "---------------"

PNFS_OUTDIR=/pnfs/dune/persistent/users/${GRID_USER}/${PNFS_PATH_APPEND}
echo "output dir is ${PNFS_OUTDIR}"

voms-proxy-info --all

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup root v6_06_08 -q e10:nu:prof
setup ifdhc v2_1_0

ups active

export IFDH_CP_UNLINK_ON_ERROR=1;
export IFDH_CP_MAXRETRIES=1;
# export IFDH_DEBUG=0;

ifdh ls ${PNFS_OUTDIR}/flux

if [ $? -ne 0 ]; then
  echo "Unable to read ${PNFS_OUTDIR}/flux. Make sure that you have created this directory and given it group write permission (chmod g+w ${PNFS_OUTDIR})."
  exit 10
fi

ifdh ls ${PNFS_OUTDIR}/logs

if [ $? -ne 0 ]; then
  echo "Unable to read ${PNFS_OUTDIR}/logs. Make sure that you have created this directory and given it group write permission (chmod g+w ${PNFS_OUTDIR})."
  exit 11
fi

mkdir inputs
cd inputs

echo "Copying inputs @ $(date)"

for FP_IT in ${INPUT_FILE_PATHS[@]}; do
  FN=${FP_IT##*/}
  echo "[INFO]: File ${FP_IT} => ${FN}"

  echo "ifdh ls ${FP_IT}"
  ifdh ls ${FP_IT}

  if [ $? -ne 0 ]; then
    echo "Unable to read input file: \"${FP_IT}\". Does it exist?"
    exit 12
  fi

  echo "[FILE: ${FP_IT}]: ifdh cp $IFDH_OPTION ${FP_IT} ${FN}"
  ifdh cp $IFDH_OPTION ${FP_IT} ${FN}

  if [ $? -ne 0 ]; then
    echo "Unable to copy input file: \"${FP_IT}\". Does it exist?"
    exit 13
  fi
done

cd ../

echo "---ls inputs---"
ls  -lah inputs
echo "---------------"

echo "------ls-------"
ls -lah
echo "---------------"

OUT_FILENAME=Fluxes.${CLUSTER}.${PROCESS}.root

echo "[INFO]: Writing output to: ${OUT_FILENAME} "

echo "Building fluxes @ $(date)"

echo "./dp_BuildFluxes -i \"inputs/*dk2nu*root\" -o ${OUT_FILENAME} -z ${DET_DIST_CM} -vb ${BINNING_DESCRIPTOR} ${RUPARG} ${FLUX_WINDOW_DESCRIPTOR} ${DK2NU_LITE_ARG} ${PPFX_ARG} &> dp_BuildFluxes.${CLUSTER}.${PROCESS}.log"
./dp_BuildFluxes -i "inputs/*dk2nu*root" -o ${OUT_FILENAME} -z ${DET_DIST_CM} -vb ${BINNING_DESCRIPTOR} ${RUPARG} ${FLUX_WINDOW_DESCRIPTOR} ${DK2NU_LITE_ARG} ${PPFX_ARG} &> dp_BuildFluxes.${CLUSTER}.${PROCESS}.log

echo "Finished."

echo "Copying output @ $(date)"

if [ ! -e ${OUT_FILENAME} ]; then
  echo "[ERROR]: Failed to produce expected output file."
  exit 1
fi

echo "ifdh cp -D $IFDH_OPTION ${OUT_FILENAME} ${PNFS_OUTDIR}/flux/"
ifdh cp -D $IFDH_OPTION ${OUT_FILENAME} ${PNFS_OUTDIR}/flux/
echo "ifdh cp -D $IFDH_OPTION dp_BuildFluxes.${CLUSTER}.${PROCESS}.log ${PNFS_OUTDIR}/logs/"
ifdh cp -D $IFDH_OPTION dp_BuildFluxes.${CLUSTER}.${PROCESS}.log ${PNFS_OUTDIR}/logs/


echo "All stop @ $(date)"
