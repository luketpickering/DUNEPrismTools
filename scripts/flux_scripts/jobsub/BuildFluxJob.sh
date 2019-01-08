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

PNFS_PATH_APPEND=""
if [ ! -z ${1} ]; then
  echo "[INFO]: Appending path to pnfs outdir: ${1}"
  PNFS_PATH_APPEND=${1}
fi

NFILES_TO_READ=1
if [ ! -z ${2} ]; then
  NFILES_TO_READ=${2}
fi

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
cp build_flux.fcl $_CONDOR_SCRATCH_DIR/

cd $_CONDOR_SCRATCH_DIR

echo "pwd is $(pwd)"
echo "------ls-------"
ls
echo "---------------"

PNFS_OUTDIR=/pnfs/dune/persistent/users/${GRID_USER}/${PNFS_PATH_APPEND}
echo "output dir is ${PNFS_OUTDIR}"

voms-proxy-info --all

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup root v6_10_08b -q e15:nu:prof
setup ifdhc v2_3_9

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

echo "------cat ./build_flux.fcl-------"
cat ./build_flux.fcl
echo "---------------"

OUT_FILENAME=Fluxes.${CLUSTER}.${PROCESS}.root

echo "[INFO]: Writing output to: ${OUT_FILENAME} "

echo "Building fluxes @ $(date)"

echo "./dp_BuildFluxes -i \"inputs/*dk2nu*root\" -o ${OUT_FILENAME} --fhicl ./build_flux.fcl &> dp_BuildFluxes.${CLUSTER}.${PROCESS}.log"
./dp_BuildFluxes -i "inputs/*dk2nu*root" -o ${OUT_FILENAME} --fhicl ./build_flux.fcl &> dp_BuildFluxes.${CLUSTER}.${PROCESS}.log

echo "Finished."

echo "------ls-------"
ls -lah
echo "---------------"

echo "Copying output @ $(date)"

echo "ifdh cp -D $IFDH_OPTION dp_BuildFluxes.${CLUSTER}.${PROCESS}.log ${PNFS_OUTDIR}/logs/"
ifdh cp -D $IFDH_OPTION dp_BuildFluxes.${CLUSTER}.${PROCESS}.log ${PNFS_OUTDIR}/logs/

if [ ! -e ${OUT_FILENAME} ]; then
  echo "[ERROR]: Failed to produce expected output file."
  exit 1
fi

echo "ifdh cp -D $IFDH_OPTION ${OUT_FILENAME} ${PNFS_OUTDIR}/flux/"
ifdh cp -D $IFDH_OPTION ${OUT_FILENAME} ${PNFS_OUTDIR}/flux/

echo "All stop @ $(date)"
