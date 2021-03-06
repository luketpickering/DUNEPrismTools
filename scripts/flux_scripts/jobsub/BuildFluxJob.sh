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

PPFX_ARGS=""
if [ ${1} != "0" ]; then
  PPFX_ARGS="--PPFX --NPPFX-Universes ${1}"
fi

if [ ${2} == "--PPFX-Components" ]; then
  PPFX_ARGS="${PPFX_ARGS} --PPFX-Components"
  shift
fi

PDGARGS=""
if [ ${2} == "--only-pdg" ]; then
  PDGARGS="--only-pdg ${3}"
  #Shift both off the end
  shift
  shift
fi

NFILES_TO_READ=1
if [ ! -z ${2} ]; then
  NFILES_TO_READ=${2}
fi

FIRST_PNFS_PATH_APPEND=""
if [ ! -z ${3} ]; then
  echo "[INFO]: Appending path to pnfs outdir: ${3}"
  FIRST_PNFS_PATH_APPEND=${3}
fi

SECOND_PNFS_PATH_APPEND=""
if [ ! -z ${4} ]; then
  echo "[INFO]: Appending path to pnfs outdir: ${4}"
  SECOND_PNFS_PATH_APPEND=${4}
fi

FIRST_FCL="ND_build_flux.fcl"
FIRST_DET="ND"
if [ ! -e ${FIRST_FCL} ]; then
  FIRST_FCL="FD_build_flux.fcl"
  FIRST_DET="FD"
else
  if [ ! -z ${SECOND_PNFS_PATH_APPEND} ]; then
    SECOND_FCL="FD_build_flux.fcl"
    SECOND_DET="FD"
  fi
fi

echo "[INFO]: Reading ${NFILES_TO_READ} files."

echo "[INFO]: JobID ${JOBSUBPARENTJOBID}, ArrayID ${JOBSUBJOBSECTION}"

(( ZERO_INDEX_LN = ${JOBSUBJOBSECTION} * ${NFILES_TO_READ} ))

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

cp dp_BuildFluxes *.so *.pcm *.h $_CONDOR_SCRATCH_DIR/
cp *build_flux.fcl $_CONDOR_SCRATCH_DIR/

cd $_CONDOR_SCRATCH_DIR

export LD_LIBRARY_PATH="$(readlink -f ./):${LD_LIBRARY_PATH}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"

echo "pwd is $(pwd)"
echo "------ls-------"
ls
echo "---------------"

voms-proxy-info --all

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup root v6_10_08b -q e15:nu:prof
setup dk2nu v01_05_01b -q e15:prof
setup ifdhc

ups active

export IFDH_CP_UNLINK_ON_ERROR=1;
export IFDH_CP_MAXRETRIES=2;
# export IFDH_DEBUG=0;

# check directories
for STAGE in FIRST SECOND; do

  GET_DET=${STAGE}_DET

  if [ -z ${!GET_DET} ]; then
    echo "[INFO]: There is no stage: ${STAGE}"
    continue
  fi

  GET_PATH_APPEND=${STAGE}_PNFS_PATH_APPEND

  if [ -z ${!GET_PATH_APPEND} ]; then
    echo "[INFO]: Expected to find output path for stage: ${STAGE}, but failed."
    exit 25
  fi

  PNFS_PATH_APPEND=${!GET_PATH_APPEND}

  PNFS_OUTDIR=/pnfs/dune/persistent/users/${GRID_USER}/${PNFS_PATH_APPEND}
  echo "Output dir for stage: ${STAGE} is ${PNFS_OUTDIR}"

  ifdh ls ${PNFS_OUTDIR}/flux &> /dev/null

  if [ $? -ne 0 ]; then
    echo "Unable to read ${PNFS_OUTDIR}/flux. Make sure that you have created this directory and given it group write permission (chmod g+w ${PNFS_OUTDIR})."
    exit 10
  fi

  ifdh ls ${PNFS_OUTDIR}/logs &> /dev/null

  if [ $? -ne 0 ]; then
    echo "Unable to read ${PNFS_OUTDIR}/logs. Make sure that you have created this directory and given it group write permission (chmod g+w ${PNFS_OUTDIR})."
    exit 11
  fi

done

echo "Writing inputs.list @ $(date)"
rm -f inputs.list; touch inputs.list

for FP_IT in ${INPUT_FILE_PATHS[@]}; do
  # echo "ifdh ls ${FP_IT}"
  # ifdh ls ${FP_IT}
  #
  # if [ $? -ne 0 ]; then
  #   echo "Unable to read input file: \"${FP_IT}\". Does it exist?"
  #   exit 12
  # fi

  XROOTFN=$( echo ${FP_IT} | sed "s|/pnfs/dune|root://fndca1.fnal.gov:1094/pnfs/fnal.gov/usr/dune|g" )
  echo "[INFO]: File ${FP_IT} => ${XROOTFN}"
  echo ${XROOTFN} >> inputs.list
done

echo "---cat inputs---"
cat inputs.list
echo "---------------"

echo "------ls-------"
ls -lah
echo "---------------"

for i in *build_flux.fcl; do
  echo "------cat ${i}-------"
  cat ${i}
  echo "---------------"
done

for STAGE in FIRST SECOND; do

  GET_DET=${STAGE}_DET

  if [ -z ${!GET_DET} ]; then
    echo "[INFO]: There is no stage: ${STAGE}"
    continue
  fi

  DET=${!GET_DET}

  FLUX_FCL=${DET}_build_flux.fcl

  OUT_FILENAME=Fluxes.${DET}.${JOBSUBPARENTJOBID}.${JOBSUBJOBSECTION}.root

  echo "[INFO]: Writing output to: ${OUT_FILENAME} "

  echo "Building fluxes @ $(date)"

  echo "./dp_BuildFluxes -i inputs.list -o ${OUT_FILENAME} --fhicl ./${FLUX_FCL} ${PPFX_ARGS} ${PDGARGS} &> dp_BuildFluxes.${JOBSUBPARENTJOBID}.${JOBSUBJOBSECTION}.log"
  ./dp_BuildFluxes -i inputs.list -o ${OUT_FILENAME} --fhicl ./${FLUX_FCL} ${PPFX_ARGS} ${PDGARGS} &> dp_BuildFluxes.${JOBSUBPARENTJOBID}.${JOBSUBJOBSECTION}.log

  echo "Finished."

  echo "------ls-------"
  ls -lah
  echo "---------------"

  GET_PATH_APPEND=${STAGE}_PNFS_PATH_APPEND
  PNFS_PATH_APPEND=${!GET_PATH_APPEND}

  PNFS_OUTDIR=/pnfs/dune/persistent/users/${GRID_USER}/${PNFS_PATH_APPEND}

  echo "Copying output @ $(date)"

  echo "ifdh cp -D $IFDH_OPTION dp_BuildFluxes.${JOBSUBPARENTJOBID}.${JOBSUBJOBSECTION}.log ${PNFS_OUTDIR}/logs/"
  ifdh cp -D $IFDH_OPTION dp_BuildFluxes.${JOBSUBPARENTJOBID}.${JOBSUBJOBSECTION}.log ${PNFS_OUTDIR}/logs/

  if [ ! -e ${OUT_FILENAME} ]; then
    echo "[ERROR]: Failed to produce expected output file."
    exit 1
  fi

  echo "ifdh cp -D $IFDH_OPTION ${OUT_FILENAME} ${PNFS_OUTDIR}/flux/"
  ifdh cp -D $IFDH_OPTION ${OUT_FILENAME} ${PNFS_OUTDIR}/flux/

done

echo "All stop @ $(date)"
