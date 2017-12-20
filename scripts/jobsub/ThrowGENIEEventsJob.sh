#!/bin/bash

if [ -z ${INPUT_TAR_FILE} ]; then
  echo "[ERROR]: Expected to recieve an input file."
  exit 1
fi

if [ ! -e inputs.list ]; then
  echo "[ERROR]: Couldn't find inputs.list"
  ls
  exit 2
fi

if [ ! -e geom.gdml ]; then
  echo "[ERROR]: Couldn't find geom.gdml"
  ls
  exit 3
fi

if [ ! -e dk2nu_FluxWindow.xml ]; then
  echo "[ERROR]: Couldn't find dk2nu_FluxWindow.xml"
  ls
  exit 4
fi

if [ -z $1 ]; then
  echo "[ERROR]: Expected to recieve a paramset name for a parameter set within dk2nu_FluxWindow.xml."
  exit 5
fi
FLUX_WINDOW_PARAMSET=${1}

DOLOGS=""
if [ ! -z ${3} ]; then
  DOLOGS="TRUE"
fi

PNFS_PATH_APPEND=""
if [ ! -z $2 ]; then
  echo "[INFO]: Appending path to pnfs outdir: $2"
  PNFS_PATH_APPEND=$2
fi

echo "JobID ${CLUSTER}, ArrayID ${PROCESS}"

(( LN = ${PROCESS} + 1 ))

INPUT_FILEPATH=$( cat inputs.list | head -${LN} | tail -1 )
INPUT_FILENAME=${INPUT_FILEPATH##*/}

echo "Working on input file: ${INPUT_FILEPATH} => ${INPUT_FILENAME}"

if [ -z ${INPUT_FILEPATH} ]; then
  echo "[ERROR]: Failed to get input file #${LN}."
  cat inputs.list
  exit 1
fi

echo "Flux window param set: ${FLUX_WINDOW_PARAMSET} "

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
  exit 6
fi

cp geom.gdml dk2nu_FluxWindow.xml $_CONDOR_SCRATCH_DIR/

cd $_CONDOR_SCRATCH_DIR

echo "pwd is $(pwd)"
echo "------ls-------"
ls
echo "---------------"

PNFS_OUTDIR=/pnfs/dune/persistent/users/${GRID_USER}/${PNFS_PATH_APPEND}
echo "output dir is ${PNFS_OUTDIR}"

voms-proxy-info --all

source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

setup larsoft v06_56_01 -q debug:e14
setup dk2nu v01_05_00j -q debug:e14
setup genie_xsec v2_12_8 -q DefaultPlusMECWithNC

setup ifdhc v2_1_0 -q debug:e14:p2713d

export IFDH_CP_UNLINK_ON_ERROR=1;
export IFDH_CP_MAXRETRIES=1;
export IFDH_DEBUG=1;

ifdh ls ${PNFS_OUTDIR}/gevgen_fnal/gntp

if [ $? -ne 0 ]; then
  echo "Unable to read ${PNFS_OUTDIR}/gevgen_fnal/gntp. Make sure that you have created this directory and given it group write permission (chmod g+w ${PNFS_OUTDIR})."
  exit 7
fi

ifdh ls ${PNFS_OUTDIR}/gevgen_fnal/rootracker

if [ $? -ne 0 ]; then
  echo "Unable to read ${PNFS_OUTDIR}/gevgen_fnal/rootracker. Make sure that you have created this directory and given it group write permission (chmod g+w ${PNFS_OUTDIR})."
  exit 7
fi

if [ ! -z ${DOLOGS} ]; then
  ifdh ls ${PNFS_OUTDIR}/gevgen_fnal/logs

  if [ $? -ne 0 ]; then
    echo "Unable to read ${PNFS_OUTDIR}/gevgen_fnal/logs. Make sure that you have created this directory and given it group write permission (chmod g+w ${PNFS_OUTDIR})."
    exit 7
  fi
fi

ifdh ls ${INPUT_FILEPATH}

if [ $? -ne 0 ]; then
  echo "Unable to read input file: \"${INPUT_FILEPATH}\". Does it exist?"
  exit 8
fi

echo "ifdh cp $IFDH_OPTION ${INPUT_FILEPATH} ${INPUT_FILENAME}"
ifdh cp $IFDH_OPTION ${INPUT_FILEPATH} ${INPUT_FILENAME}

echo "------ls-------"
ls
echo "---------------"

GNTPOUT_FILENAME=${INPUT_FILENAME%%\.dk2nu\.root}.gntpc.${CLUSTER}.${PROCESS}.root
RTOUT_FILENAME=${INPUT_FILENAME%%\.dk2nu\.root}.rootracker.${CLUSTER}.${PROCESS}.root

echo "[INFO]: Writing gntpc output to: ${GNTPOUT_FILENAME}"
echo "[INFO]: Writing rootracker output to: ${RTOUT_FILENAME}"

GXMLPATH_OLD=${GXMLPATH}

export GXMLPATH=$(pwd):${GXMLPATH}
export GDK2NUFLUXXML=dk2nu_FluxWindow.xml

echo "[INFO]: OLDGXMLPATH: ${GXMLPATH_OLD}, GXMLPATH=${GXMLPATH}."

echo "Running GENIE @ $(date)"
echo "gevgen_fnal --message-thresholds ${GENIE}/config/Messenger_whisper.xml -r ${PROCESS} -f ${INPUT_FILENAME},${FLUX_WINDOW_PARAMSET} -g geom.gdml -e 2E16 -L cm -D g_cm3 --cross-sections ${GXMLPATH_OLD}/gxspl-FNALsmall.xml &> gevgen_fnal.${CLUSTER}.${PROCESS}.log"
gevgen_fnal --message-thresholds ${GENIE}/config/Messenger_whisper.xml -r ${PROCESS} -f ${INPUT_FILENAME},${FLUX_WINDOW_PARAMSET} -g geom.gdml -e 2E16 -L cm -D g_cm3 --cross-sections ${GXMLPATH_OLD}/gxspl-FNALsmall.xml &> gevgen_fnal.${CLUSTER}.${PROCESS}.log
echo "Finished."
echo "Converting to rootracker @ $(date)"
echo "------ls-------"
ls
echo "---------------"
mv gntp.${PROCESS}.ghep.root ${GNTPOUT_FILENAME}
gntpc -i ${GNTPOUT_FILENAME} -f rootracker -o ${RTOUT_FILENAME} &> /dev/null
echo "Finished."

echo "------ls-------"
ls
echo "---------------"

echo "Copying output @ $(date)"
echo "ifdh cp -D $IFDH_OPTION ${GNTPOUT_FILENAME} ${PNFS_OUTDIR}/gevgen_fnal/gntp"
ifdh cp -D $IFDH_OPTION ${GNTPOUT_FILENAME} ${PNFS_OUTDIR}/gevgen_fnal/gntp
echo "ifdh cp -D $IFDH_OPTION ${RTOUT_FILENAME} ${PNFS_OUTDIR}/gevgen_fnal/rootracker"
ifdh cp -D $IFDH_OPTION ${RTOUT_FILENAME} ${PNFS_OUTDIR}/gevgen_fnal/rootracker
if [ ! -z ${DOLOGS} ]; then
  echo "ifdh cp -D $IFDH_OPTION gevgen_fnal.${CLUSTER}.${PROCESS}.log ${PNFS_OUTDIR}/gevgen_fnal/logs"
  ifdh cp -D $IFDH_OPTION gevgen_fnal.${CLUSTER}.${PROCESS}.log ${PNFS_OUTDIR}/gevgen_fnal/logs
fi
echo "All stop"
