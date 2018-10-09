#!/bin/sh

PNFS_PATH_APPEND=""
DK2NU_INPUT_DIR=""
NMAXJOBS=""
NPERJOB="10"
INPUTLIST=""
LIFETIME_EXP="30m"
DISK_EXP="1GB"
MEM_EXP="1GB"
INCLUDE_PPFX="0"

while [[ ${#} -gt 0 ]]; do

  key="$1"
  case $key in

      -p|--pnfs-path-append)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      PNFS_PATH_APPEND="$2"
      echo "[OPT]: Writing output to /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}"
      shift # past argument
      ;;

      -i|--dk2nu-input-directory)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      DK2NU_INPUT_DIR="$2"
      echo "[OPT]: Looking for input files in \"${DK2NU_INPUT_DIR}\"."
      shift # past argument
      ;;

      -I|--dk2nu-input-file-list)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      INPUTLIST=$(readlink -f $2)
      echo "[OPT]: Using \"${INPUTLIST}\" as a dk2nu input file list."
      shift # past argument
      ;;

      -P|--PPFX)

      INCLUDE_PPFX="1"
      echo "[OPT]: Looking for PPFX friend files."
      ;;

      -n|--n-per-job)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      NPERJOB="$2"
      echo "[OPT]: Each job will handle \"${NPERJOB}\" input files."
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

      --expected-walltime)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      LIFETIME_EXP="$2"
      echo "[OPT]: Expecting a run time of \"${LIFETIME_EXP}\"."
      shift # past argument
      ;;

      --expected-disk)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      DISK_EXP="$2"
      echo "[OPT]: Expecting to use \"${DISK_EXP}\" node disk space."
      shift # past argument
      ;;


      --expected-mem)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      MEM_EXP="$2"
      echo "[OPT]: Expecting a maximum of \"${MEM_EXP}\" memory usage."
      shift # past argument
      ;;

      -?|--help)

      echo "[RUNLIKE] ${SCRIPTNAME}"
      echo -e "\t-p|--pnfs-path-append      : Path to append to output path: /pnfs/dune/persistent/users/${USER}/"
      echo -e "\t-d|--detector-distance     : Detector distance from beam z=0 in cm."
      echo -e "\t-i|--dk2nu-input-directory : Input directory to search for dk2nu.root files"
      echo -e "\t-I|--dk2nu-input-file-list : Newline separated list of files to use as input. Must be located on dcache."
      echo -e "\t-P|--PPFX                  : Look for PPFX friend trees to input dk2nu files."
      echo -e "\t-n|--n-per-job             : Number of files to run per job. (default: 10)"
      echo -e "\t-N|--NMAXJobs              : Maximum number of jobs to submit."
      echo -e "\t--expected-disk            : Expected disk usage to pass to jobsub -- approx 100MB* the value passed to -\'n\' (default: 1GB)"
      echo -e "\t--expected-mem             : Expected mem usage to pass to jobsub (default: 2GB)"
      echo -e "\t--expected-walltime        : Expected disk usage to pass to jobsub (default: 20m)"
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

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh

setup jobsub_client
setup ifdhc

if [ -e sub_tmp ]; then rm -r sub_tmp; fi

mkdir sub_tmp
cd sub_tmp

# export IFDH_DEBUG=1

if [ ! "${EXPERIMENT}" ]; then
  echo "[ERROR]: No experiment. Expected \${EXPERIMENT} to be defined."
  exit 1
fi

if [ "${EXPERIMENT}" != "dune" ]; then
  echo "[WARN]: \${EXPERIMENT} != \"dune\", \${EXPERIMENT} = \"${EXPERIMENT}\"."
fi

if [ ! "${DUNEPRISMTOOLSROOT}" ]; then
  echo "[ERROR]: \$DUNEPRISMTOOLSROOT was not defined, you must set up DUNEPrismTools before deploying."
  exit 1
fi

NINPUTS=0
if [ ! -z ${INPUTLIST} ]; then
  if [ ! -e "${INPUTLIST}" ]; then
    echo "[ERROR]: No or non-existant dk2nu input list: \"${INPUTLIST}\"."
    exit 1
  fi
  cp ${INPUTLIST} inputs.list

  NINPUTS=$(cat inputs.list | wc -l)

  if [ -z ${NINPUTS} ]; then
    echo "[ERROR]: Found no inputs in: ${INPUTLIST}."
    exit 1
  fi

else

  if [ ! "${DK2NU_INPUT_DIR}" ] || [ ! -e "${DK2NU_INPUT_DIR}" ]; then
    echo "[ERROR]: No or non-existant dk2nu input dir: \"${DK2NU_INPUT_DIR}\"."
    exit 1
  fi

  PNFS_PREFIX=${DK2NU_INPUT_DIR%%/dune*}

  if [ "${PNFS_PREFIX}" == "/pnfs" ]; then
    ifdh ls ${DK2NU_INPUT_DIR} | grep "dk2nu.root$" > inputs.list
  else
    echo -e "import os,re\nfor x in os.listdir('${DK2NU_INPUT_DIR}'):\n if re.match('.*dk2nu.root',x):\n  print x;" | python > inputs.list
  fi

  NINPUTS=$(cat inputs.list | wc -l)

  if [ -z ${NINPUTS} ]; then
    echo "[ERROR]: Found no inputs in: ${DK2NU_INPUT_DIR}."
    exit 1
  fi
fi

if [ "${INCLUDE_PPFX}" == "1"  ]; then

  rm -f ppfx_inputs.list
  touch ppfx_inputs.list

  for INP in cat inputs.list; do
    IDIR=${INP%/*}
    PPFX_IDIR=$(echo ${IDIR} | sed 's:neutrino/flux:neutrino/ppfx/hadron/flux:g')
    INPF=${INP##*/}
    PPFX_INPF=$(echo ${INPF} | sed 's:neutrino_:neutrino_ppfx_friend_:g' | sed 's:dk2nu\.::g')
    echo ${PPFX_IDIR}/${PPFX_INPF} >> ppfx_inputs.list
  done

fi

NJOBSTORUN=$(python -c "from math import ceil; print int(ceil(float(${NINPUTS})/float(${NPERJOB}))); ")
if [ ! -z ${NMAXJOBS} ]; then
  NJOBSTORUN=$(python -c "print min(${NMAXJOBS},${NJOBSTORUN})")
fi

echo "[INFO]: Found ${NINPUTS} inputs, will run ${NJOBSTORUN} jobs."

if [ ! -e ${DUNEPRISMTOOLSROOT}/bin/dp_MakeLitedk2nu ]; then
  echo "[ERROR]: It appears that DUNEPrismTools was not built, expected to find \"${DUNEPRISMTOOLSROOT}/bin/dp_MakeLitedk2nu\"."
  exit 1
fi

ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/slimflux

if [ $? -ne 0 ]; then
  echo "Attempting to make /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/slimflux."
  mkdir -p /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/slimflux
  ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/slimflux
  if [ $? -ne 0 ]; then
    echo "Unable to make /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/slimflux."
    exit 7
  fi
  echo "Made /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/slimflux."
fi

cp ${DUNEPRISMTOOLSROOT}/bin/dp_MakeLitedk2nu .

tar -zcvf apps.@DUNEPrismTools_VERSION_STRING@.tar.gz dp_MakeLitedk2nu *inputs.list

echo "Submitting job: jobsub_submit --group=${EXPERIMENT} --jobid-output-only --resource-provides=usage_model=OPPORTUNISTIC --expected-lifetime=${LIFETIME_EXP} --disk=${DISK_EXP} -N ${NJOBSTORUN} --memory=${MEM_EXP} --cpu=1 --OS=SL6 --tar_file_name=dropbox://apps.@DUNEPrismTools_VERSION_STRING@.tar.gz file://${DUNEPRISMTOOLSROOT}/scripts/SlimFluxJob.sh ${PNFS_PATH_APPEND} ${NPERJOB}"
JID=$(jobsub_submit --group=${EXPERIMENT} --jobid-output-only --resource-provides=usage_model=OPPORTUNISTIC --expected-lifetime=${LIFETIME_EXP} --disk=${DISK_EXP} -N ${NJOBSTORUN} --memory=${MEM_EXP} --cpu=1 --OS=SL6 --tar_file_name=dropbox://apps.@DUNEPrismTools_VERSION_STRING@.tar.gz file://${DUNEPRISMTOOLSROOT}/scripts/SlimFluxJob.sh ${PNFS_PATH_APPEND} ${NPERJOB})
echo "Done."

cd ../
rm -r sub_tmp

echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; jobsub_q --jobid=${JID}" > JID_${JID}.q.sh
echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; jobsub_rm --jobid=${JID}" > JID_${JID}.rm.sh
echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; mkdir ${JID}_log; cd ${JID}_log; jobsub_fetchlog --jobid=${JID}; tar -zxvf *.tgz" > JID_${JID}.fetchlog.sh

chmod +x JID_${JID}.q.sh JID_${JID}.rm.sh JID_${JID}.fetchlog.sh
