#!/bin/bash --login
#PBS -l walltime=1:59:00
#PBS -l mem=512MB
#PBS -j oe

touch ${PBS_O_WORKDIR}/CondenseProcess.${PBS_ARRAYID}.running

if [ -z ${ENVSETUPSCRIPT} ]; then
  echo "[ERROR]: Was not passed the location of an environment setup script."
  exit 1
fi

source ${ENVSETUPSCRIPT}

echo Job $PBS_JOBNAME started on $HOST at $(date "+%Y.%m.%d %H:%M:%S %Z")

if [ -z ${INPUT_FILE_LIST} ]; then
  echo "[ERROR]: Was not passed a list of files to work on."
  exit 1
fi

if [ -z ${RUNPLAN_CONFIG} ]; then
  echo "[ERROR]: Was not passed a run plan."
  exit 1
fi

if [ -z ${CONDENSED_OUTPUT_DIR} ]; then
  echo "[ERROR]: Was not passed a condensed file output directory."
  exit 1
fi

if [ -z ${PROCESSED_OUTPUT_DIR} ]; then
  echo "[ERROR]: Was not passed a processed file output directory."
  exit 1
fi

if [ -z ${CONDENSERCONFIG} ]; then
  echo "[ERROR]: Was not passed a condenser configuration."
  exit 1
fi

PPFARG=""
if [ ! -z ${POTPERFILE} ]; then
  echo "[INFO]: Using ${POTPERFILE} as POT per file."
  PPFARG=" -P ${POTPERFILE} "
fi

UENCCC=$(python -c "import urllib; print urllib.unquote('''${CONDENSERCONFIG}''')")

echo "[INFO]: Unencoded CONDENSERCONFIG = \"${UENCCC}\""

G4PYFile=$(cat ${INPUT_FILE_LIST} | head -${PBS_ARRAYID} | tail -1 | cut -d " " -f 1)
RTFile=$(cat ${INPUT_FILE_LIST} | head -${PBS_ARRAYID} | tail -1 | cut -d " " -f 2)

mkdir ${TMPDIR}/CondenseAndProcess_${PBS_JOBID}.${PBS_ARRAYID}

cd ${TMPDIR}/CondenseAndProcess_${PBS_JOBID}.${PBS_ARRAYID}

cp ${G4PYFile} ${RTFile} ${RUNPLAN_CONFIG} ./

G4PYFileName=${G4PYFile##*/}
RTFileName=${RTFile##*/}
RPFileName=${RUNPLAN_CONFIG##*/}

CONDFileName=$( echo ${G4PYFileName} | awk '{split($0,a,"."); print a[1] "." a[2] ".Condensed.root"};' )

echo "Condensing at $(date "+%Y.%m.%d %H:%M:%S %Z")"
echo "dp_Condenser -i ${G4PYFileName} -ir ${RTFileName} ${UENCCC} -o ${CONDFileName} ${PPFARG}"
dp_Condenser -i ${G4PYFileName} -ir ${RTFileName} ${UENCCC} -o ${CONDFileName} ${PPFARG}

ProcFileName=$( echo ${G4PYFileName} | awk '{split($0,a,"."); print a[1] "." a[2] ".Processed.root"};' )

echo "Processing at $(date "+%Y.%m.%d %H:%M:%S %Z")"
echo "dp_StopProcessor -i ${CONDFileName} -r ${RPFileName} -o ${ProcFileName}"
dp_StopProcessor -i ${CONDFileName} -r ${RPFileName} -o ${ProcFileName}

echo "Copy back at $(date "+%Y.%m.%d %H:%M:%S %Z")"
echo "cp ${CONDFileName} ${CONDENSED_OUTPUT_DIR}"
cp ${CONDFileName} ${CONDENSED_OUTPUT_DIR}

echo "cp ${ProcFileName} ${PROCESSED_OUTPUT_DIR}"
cp ${ProcFileName} ${PROCESSED_OUTPUT_DIR}

echo "Job finished at $(date "+%Y.%m.%d %H:%M:%S %Z")"

rm ${PBS_O_WORKDIR}/CondenseProcess.${PBS_ARRAYID}.running
