#!/bin/bash

COPYBACKINTERMEDIATES=0

if [ -z ${1} ]; then
  echo "[ERROR]: Expected to recieve a macro name to run."
  exit 7
fi
BEAM_MACRO=${1}

if [ ! -e ${BEAM_MACRO} ]; then
  echo "[ERROR]: Couldn't find G4 macro: \"${BEAM_MACRO}\""
  ls
  exit 8
fi

PNFS_PATH_APPEND=""
if [ ! -z ${2} ]; then
  echo "[INFO]: Appending path to pnfs outdir: ${2}"
  PNFS_PATH_APPEND=${2}
fi

USE_PPFX="0"
if [ ! -z ${3} ]; then
  USE_PPFX=${3}
fi
echo "[INFO]: Using PPFX ? ${USE_PPFX}"

USE_ALL_PPFX_WEIGHTS="0"
if [ ! -z ${4} ]; then
  USE_ALL_PPFX_WEIGHTS=${4}
fi
echo "[INFO]: Using all PPFX weights? ${USE_ALL_PPFX_WEIGHTS}"


if [ -z ${INPUT_TAR_FILE} ]; then
  echo "[ERROR]: Expected to recieve an input file."
  exit 1
fi

APPS="g4lbnf dp_MakeLitedk2nu"
DPLIBS="libTH2Jagged.so"
if [ "${USE_PPFX}" == "1" ]; then
  APPS="g4lbnf Make_dk2nu_FriendTree dp_MakeLitedk2nu"
fi

for app in ${APPS}; do
  if [ ! -e ${app} ]; then
    echo "[ERROR]: Couldn't find expected ${app}"
    exit 2
  fi
done
for lib in ${DPLIBS}; do
  if [ ! -e ${lib} ]; then
    echo "[ERROR]: Couldn't find expected DUNEPrismTools dependency: ${lib}"
    exit 3
  fi
done

G4LBNE_LIBS="libmyVersion.so libg4lbnfDict.so g4lbnfCint_rdict.pcm"
G4LBNE_INPUTS="gdml locations"
for lib in ${G4LBNE_LIBS}; do
  if [ ! -e ${lib} ]; then
    echo "[ERROR]: Couldn't find expected g4lbnf dependency: ${lib}"
    exit 3
  fi
done
for dir in ${G4LBNE_INPUTS}; do
  if [ ! -e ${dir} ]; then
    echo "[ERROR]: Couldn't find expected g4lbnf inputs directory: ${dir}"
    exit 4
  fi
done

PPFX_LIBS=""
PPFX_INPUTS=""
if [ "${USE_PPFX}" == "1" ]; then
  PPFX_LIBS="libppfx.so"
  for lib in ${PPFX_LIBS}; do
    if [ ! -e ${lib} ]; then
      echo "[ERROR]: Couldn't find expected ppfx dependency: ${lib}"
      exit 4
    fi
  done

  PPFX_INPUTS="ppfx/data ppfx/uncertainties"
  for dir in ${PPFX_INPUTS}; do
    if [ ! -e ${dir} ]; then
      echo "[ERROR]: Couldn't find expected ppfx inputs directory: ${dir}"
      exit 4
    fi
  done

  if [ ! -e ppfx.xml ]; then
    echo "[ERROR]: Couldn't find expected ppfx inputs file: ppfx.xml"
    exit 6
  fi
fi

if [ ! -e setups/setup.sh ]; then
  echo "[ERROR]: Couldn't find expected setup script: setups/setup.sh"
  exit 5
fi


echo "[INFO]: JobID ${CLUSTER}, ArrayID ${PROCESS}"
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

mv ${APPS} ${DPLIBS} ${G4LBNE_LIBS} ${G4LBNE_INPUTS} setups ${BEAM_MACRO} $_CONDOR_SCRATCH_DIR/

if [ "${USE_PPFX}" == "1" ]; then
  mv ${PPFX_LIBS} ppfx.xml $_CONDOR_SCRATCH_DIR/
  mv ppfx $_CONDOR_SCRATCH_DIR/
fi

cd $_CONDOR_SCRATCH_DIR

mkdir lib
mv *.so lib/
export LD_LIBRARY_PATH="$(readlink -f lib):${LD_LIBRARY_PATH}"
echo "LD_LIBRARY_PATH=${LD_LIBRARY_PATH}"

echo "pwd is $(pwd)"
echo "------ls-------"
ls -lah
echo "---------------"

PNFS_OUTDIR=/pnfs/dune/persistent/users/${GRID_USER}/${PNFS_PATH_APPEND}
echo "Output dir is ${PNFS_OUTDIR}"

voms-proxy-info --all

source setups/setup.sh

echo "${G4LIB}"
echo "------ls-------"
ls ${G4LIB}
echo "---------------"


export IFDH_CP_UNLINK_ON_ERROR=1;
export IFDH_CP_MAXRETRIES=1;

ifdh ls ${PNFS_OUTDIR}/dk2nulite

if [ $? -ne 0 ]; then
  echo "Unable to read ${PNFS_OUTDIR}/dk2nulite. Make sure that you have created this directory and given it group write permission (chmod g+w ${PNFS_OUTDIR})."
  exit 10
fi

ifdh ls ${PNFS_OUTDIR}/logs

if [ $? -ne 0 ]; then
  echo "Unable to read ${PNFS_OUTDIR}/logs. Make sure that you have created this directory and given it group write permission (chmod g+w ${PNFS_OUTDIR})."
  exit 11
fi

G4LBNE_DK2NUFILE_PREFIX=g4lbnf

cp ${BEAM_MACRO} ${BEAM_MACRO}.inp

SEED=$(od -A n -t dI -N 2 /dev/urandom | tr -d " ")

cat ${BEAM_MACRO}.inp | sed "s/__RUNID__/${PROCESS}/g" | sed "s/__SEED__/${SEED}/g" | sed "s/__DK2NUFILENAME__/${G4LBNE_DK2NUFILE_PREFIX}/g" > ${BEAM_MACRO}

G4LBNE_DK2NUFILE=${G4LBNE_DK2NUFILE_PREFIX}_$(printf "%05g" ${PROCESS}).dk2nu.root

echo "====START====MACRO============"
cat ${BEAM_MACRO}
echo "=====END=====MACRO============"

echo "./g4lbnf --physicslist QGSP_BERT ${BEAM_MACRO} &> g4lbnf.${CLUSTER}.${PROCESS}.log @ $(date)"
./g4lbnf --physicslist QGSP_BERT ${BEAM_MACRO} &> g4lbnf.${CLUSTER}.${PROCESS}.log
echo "Finished @ $(date)"

find . -name "LBNElocations.txt"
echo $(pwd)

echo "pwd is $(pwd)"
echo "------ls-------"
ls -lah
echo "---------------"

if [ -e  g4lbnf.${CLUSTER}.${PROCESS}.log ]; then
  echo "ifdh cp -D $IFDH_OPTION  g4lbnf.${CLUSTER}.${PROCESS}.log ${PNFS_OUTDIR}/logs/"
  ifdh cp -D $IFDH_OPTION  g4lbnf.${CLUSTER}.${PROCESS}.log ${PNFS_OUTDIR}/logs/
fi

for i in $(ls core.*); do
  echo "[ERROR]: Found core file: ${i}."
  echo "ifdh cp -D $IFDH_OPTION core.* ${PNFS_OUTDIR}/logs/"
  ifdh cp -D $IFDH_OPTION core.* ${PNFS_OUTDIR}/logs/
  exit 11
done

if [ ! -e ${G4LBNE_DK2NUFILE} ]; then
  echo "g4lbnf failed to produce a dk2nu file. Cannot continue."
  exit 12
fi

if [ "${COPYBACKINTERMEDIATES}" == "1" ]; then
    echo "ifdh cp -D $IFDH_OPTION  ${G4LBNE_DK2NUFILE} ${PNFS_OUTDIR}/dk2nulite/"
    ifdh cp -D $IFDH_OPTION  ${G4LBNE_DK2NUFILE} ${PNFS_OUTDIR}/dk2nulite/
fi

PPFX_FRIENDFILE=""
if [ "${USE_PPFX}" == "1" ]; then
  PPFX_FRIENDFILE=ppfx.dk2nu_friend.${CLUSTER}.${PROCESS}.root

  export MODE="OPT"
  export PPFX_DIR=$(readlink -f ppfx)

  echo "------ldd-------"
  ldd Make_dk2nu_FriendTree
  echo "---------------"

  ALLWEIGHTSCMD=""
  if [ "${USE_ALL_PPFX_WEIGHTS}" == "1" ]; then
    ALLWEIGHTSCMD="AllWeights"
  fi

  echo "./Make_dk2nu_FriendTree ${G4LBNE_DK2NUFILE} ${PPFX_FRIENDFILE} ppfx.xml ${ALLWEIGHTSCMD} &> ppfx.${CLUSTER}.${PROCESS}.log @ $(date)"
  ./Make_dk2nu_FriendTree ${G4LBNE_DK2NUFILE} ${PPFX_FRIENDFILE} ppfx.xml ${ALLWEIGHTSCMD} &> ppfx.${CLUSTER}.${PROCESS}.log
  echo "Finished @ $(date)"

  echo "pwd is $(pwd)"
  echo "------ls-------"
  ls
  echo "---------------"

  if [ -e ppfx.${CLUSTER}.${PROCESS}.log ]; then
    echo "ifdh cp -D $IFDH_OPTION ppfx.${CLUSTER}.${PROCESS}.log ${PNFS_OUTDIR}/logs/"
    ifdh cp -D $IFDH_OPTION ppfx.${CLUSTER}.${PROCESS}.log ${PNFS_OUTDIR}/logs/
  fi

  for i in $(ls core.*); do
    echo "[ERROR]: Found core file: ${i}."
    echo "ifdh cp -D $IFDH_OPTION core.* ${PNFS_OUTDIR}/logs/"
    ifdh cp -D $IFDH_OPTION core.* ${PNFS_OUTDIR}/logs/
    exit 11
  done

  if [ ! -e ${PPFX_FRIENDFILE} ]; then
    echo "ppfx failed to produce a dk2nu friend file. Cannot continue."
    exit 13
  fi

  if [ "${COPYBACKINTERMEDIATES}" == "1" ]; then
      echo "ifdh cp -D $IFDH_OPTION ${PPFX_FRIENDFILE} ${PNFS_OUTDIR}/dk2nulite/"
      ifdh cp -D $IFDH_OPTION ${PPFX_FRIENDFILE} ${PNFS_OUTDIR}/dk2nulite/
  fi
fi

DK2NULITEFILE=g4lbnf.dk2nulite.${CLUSTER}.${PROCESS}.root

if [ "${USE_PPFX}" == "1" ]; then
  echo "./dp_MakeLitedk2nu -i ${G4LBNE_DK2NUFILE} -p ${PPFX_FRIENDFILE} -o ${DK2NULITEFILE}  @ $(date)"
  ./dp_MakeLitedk2nu -i ${G4LBNE_DK2NUFILE} -p ${PPFX_FRIENDFILE} -o ${DK2NULITEFILE}
else
  echo "./dp_MakeLitedk2nu -i ${G4LBNE_DK2NUFILE} -o ${DK2NULITEFILE}  @ $(date)"
  ./dp_MakeLitedk2nu -i ${G4LBNE_DK2NUFILE} -o ${DK2NULITEFILE}
fi

echo "Finished @ $(date)"

echo "pwd is $(pwd)"
echo "------ls-------"
ls -lah
echo "---------------"

if [ ! -e ${DK2NULITEFILE} ]; then
  echo "dp_MakeLitedk2nu failed to produce a dk2nu lite file. Job failed."
  exit 13
else
  echo "ifdh cp -D $IFDH_OPTION ${DK2NULITEFILE} ${PNFS_OUTDIR}/dk2nulite/"
  ifdh cp -D $IFDH_OPTION ${DK2NULITEFILE} ${PNFS_OUTDIR}/dk2nulite/
fi

echo "All stop @ $(date)"
