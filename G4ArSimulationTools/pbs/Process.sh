#!/bin/bash --login
#PBS -l walltime=0:20:00
#PBS -l mem=512MB
#PBS -j oe

touch ${PBS_O_WORKDIR}/Process.${PBS_ARRAYID}.running

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

if [ -z ${PROCESSED_OUTPUT_DIR} ]; then
  echo "[ERROR]: Was not passed a processed file output directory."
  exit 1
fi

CONDFILE=$(cat ${INPUT_FILE_LIST} | head -${PBS_ARRAYID} | tail -1)

mkdir ${TMPDIR}/Process_${PBS_JOBID}.${PBS_ARRAYID}

cd ${TMPDIR}/Process_${PBS_JOBID}.${PBS_ARRAYID}

cp ${CONDFILE} ${RUNPLAN_CONFIG} ./

CONDFileName=${CONDFILE##*/}
RPFileName=${RUNPLAN_CONFIG##*/}

ProcFileName=$( echo ${CONDFileName} | sed "s/Condensed/Processed/g" )

echo "Processing at $(date "+%Y.%m.%d %H:%M:%S %Z")"
echo "dp_FullDetTreeStopProcessor -i ${CONDFileName} -r ${RPFileName} -o ${ProcFileName}"
dp_FullDetTreeStopProcessor -i ${CONDFileName} -r ${RPFileName} -o ${ProcFileName}

echo "cp ${ProcFileName} ${PROCESSED_OUTPUT_DIR}"
cp ${ProcFileName} ${PROCESSED_OUTPUT_DIR}

echo "Job finished at $(date "+%Y.%m.%d %H:%M:%S %Z")"

rm ${PBS_O_WORKDIR}/Process.${PBS_ARRAYID}.running
