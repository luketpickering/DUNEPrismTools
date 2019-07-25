#!/bin/sh

PNFS_PATH_APPEND=""
DK2NU_INPUT_DIR=""
NMAXJOBS=""
NPERJOB="10"
NFILESKIP="0"
INPUTLIST=""
LIFETIME_EXP="30m"
DISK_EXP="1GB"
MEM_EXP="2GB"
FORCE_REMOVE="0"
CONFIG_FCL_ND="DUMMY.fcl"
CONFIG_FCL_FD="DUMMY.fcl"
DO_ND="0"
DO_FD="0"
ONLY_PDG="0"

NPPFXU="0"
USE_PPFX_COMPONENTS=""

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

      -n|--n-per-job)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      NPERJOB="$2"
      echo "[OPT]: Each job will handle \"${NPERJOB}\" input files."
      shift # past argument
      ;;

      -FN|--fhicl-ND)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      CONFIG_FCL_ND="$2"
      DO_ND="1"
      echo "[OPT]: Using ${CONFIG_FCL_ND} for ND fhicl."
      shift # past argument
      ;;

      -FF|--fhicl-FD)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      CONFIG_FCL_FD="$2"
      DO_FD="1"
      echo "[OPT]: Using ${CONFIG_FCL_FD} for FD fhicl."
      shift # past argument
      ;;

      -s|--n-files-skip)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      NFILESKIP="$2"
      echo "[OPT]: Will ignore the first \"${NFILESKIP}\" input files."
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

      --NPPFXU)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      NPPFXU="$2"
      echo "[OPT]: Expecting \"${NPPFXU}\" PPFX Universes."
      shift # past argument
      ;;

      --PPFX-Components)
      USE_PPFX_COMPONENTS="--PPFX-Components"
      echo "[OPT]: Building PPFX component tweak predictions"
      ;;

      --only-pdg)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      ONLY_PDG="$2"
      echo "[OPT]: Only building predictions for nu_pdg=\"${ONLY_PDG}\"."
      shift # past argument
      ;;

      -f|--force-remove)

      FORCE_REMOVE="1"
      echo "[OPT]: Will remove output directories if they exist."
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
      echo -e "\t-i|--dk2nu-input-directory : Input directory to search for dk2nu.root files"
      echo -e "\t-I|--dk2nu-input-file-list : Newline separated list of files to use as input. Must be located on dcache."
      echo -e "\t-n|--n-per-job             : Number of files to run per job. (default: 10)"
      echo -e "\t-s|--n-files-skip          : Number of files to to skip before building jobs. (default: 10)"
      echo -e "\t-FN|--fhicl-ND             : Flux calculation configuration file for ND."
      echo -e "\t-FD|--fhicl-FD             : Flux calculation configuration file for FD."
      echo -e "\t--NPPFXU <#ppfx universes> : Number of PPFX universes to use."
      echo -e "\t--PPFX-Components          : Build predictions for each PPFX universe for each reweight group."
      echo -e "\t--only-pdg                 : Only build predictions for specified neutrino species."
      echo -e "\t-f|--force-remove          : If output directories already exist, force remove them."
      echo -e "\t--expected-disk            : Expected disk usage to pass to jobsub -- approx 100MB* the value passed to -\'n\' (default: 1GB)"
      echo -e "\t--expected-mem             : Expected mem usage to pass to jobsub -- Scales with the number of detector stops in the xml passed to \'--r\' (default: 2GB)"
      echo -e "\t--expected-walltime        : Expected disk usage to pass to jobsub -- Scales with both of the above, but 20m for 1 million entries and 1000 fluxes to build is about right (default: 20m)"
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

  if [ ! -z ${NFILESKIP} ]; then
    NINPUTS=$(cat inputs.list | wc -l)
    N=$((NINPUTS - NFILESKIP))
    mv inputs.list inputs.list.full
    cat inputs.list.full | tail -n ${N} > inputs.list
  fi

  NINPUTS=$(cat inputs.list | wc -l)

  if [ "${NINPUTS}" == "0" ]; then
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
    voms-proxy-info --all

    echo "[INFO]: ifdh ls ${DK2NU_INPUT_DIR}"
    ifdh ls ${DK2NU_INPUT_DIR} | grep "dk2nu.*root$" > inputs.list
  else
    echo -e "import os,re\nfor x in os.listdir('${DK2NU_INPUT_DIR}'):\n if re.match('.*dk2nu.*root',x):\n  print x;" | python > inputs.list
  fi

  if [ ! -z ${NFILESKIP} ]; then
    NINPUTS=$(cat inputs.list | wc -l)
    N=$((NINPUTS - NFILESKIP))
    mv inputs.list inputs.list.full
    cat inputs.list.full | tail -n ${N} > inputs.list
  fi

  NINPUTS=$(cat inputs.list | wc -l)

  if [ "${NINPUTS}" == "0" ]; then
    echo "[ERROR]: Found no inputs in: ${DK2NU_INPUT_DIR}."
    exit 1
  fi
fi

NJOBSTORUN=$(python -c "from math import ceil; print int(ceil(float(${NINPUTS})/float(${NPERJOB}))); ")
if [ ! -z ${NMAXJOBS} ]; then
  NJOBSTORUN=$(python -c "print min(${NMAXJOBS},${NJOBSTORUN})")
fi

echo "[INFO]: Found ${NINPUTS} inputs, will run ${NJOBSTORUN} jobs."

if [ ! -e ${DUNEPRISMTOOLSROOT}/bin/dp_BuildFluxes ]; then
  echo "[ERROR]: It appears that DUNEPrismTools was not built, expected to find \"${DUNEPRISMTOOLSROOT}/bin/dp_BuildFluxes\"."
  exit 1
fi

if [ -e ${DUNEPRISMTOOLSROOT}/fcl/${CONFIG_FCL_ND} ]; then
  ND_PATH_APPEND=$(echo ${PNFS_PATH_APPEND} | sed "s/__DET__/ND/g")
  echo "[INFO]: Doing ND."
fi

if [ -e ${DUNEPRISMTOOLSROOT}/fcl/${CONFIG_FCL_FD} ]; then
  FD_PATH_APPEND=$(echo ${PNFS_PATH_APPEND} | sed "s/__DET__/FD/g")
  echo "[INFO]: Doing FD."
fi

for DET in ND FD; do
  CHECK_VAR=DO_${DET}

  if [ ${!CHECK_VAR} != "1" ]; then
    echo "[INFO]: Not processing ${DET} this run."
    continue
  fi

  GET_PATH_APPEND=${DET}_PATH_APPEND
  DET_PATH_APPEND=${!GET_PATH_APPEND}

  ifdh ls /pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/flux &> /dev/null

  if [ $? -ne 0 ]; then
    mkdir -p /pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/flux
    ifdh ls /pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/flux &> /dev/null
    if [ $? -ne 0 ]; then
      echo "Unable to make /pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/flux."
      exit 7
    fi
  elif [ ${FORCE_REMOVE} == "1" ]; then
    echo "[INFO]: Force removing previous existant output directories: \"/pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/flux\" "
    rm -rf /pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/flux
    mkdir -p /pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/flux
  fi

  ifdh ls /pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/logs &> /dev/null

  if [ $? -ne 0 ]; then
    mkdir -p /pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/logs
    ifdh ls /pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/logs &> /dev/null
    if [ $? -ne 0 ]; then
      echo "Unable to make /pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/logs."
      exit 7
    fi
  elif [ ${FORCE_REMOVE} == "1" ]; then
    echo "[INFO]: Force removing previous existant output directories: \"/pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/logs\" "
    rm -rf /pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/logs
    mkdir -p /pnfs/dune/persistent/users/${USER}/${DET_PATH_APPEND}/logs
  fi
done # End loop over det

cp ${DUNEPRISMTOOLSROOT}/bin/dp_BuildFluxes .
cp ${DUNEPRISMTOOLSROOT}/lib/libTH2Jagged.so ./

if [ ${DO_ND} == "1" ]; then
  if [ ! -e ${DUNEPRISMTOOLSROOT}/fcl/${CONFIG_FCL_ND} ]; then
    echo "[ERROR]: Expected to find ${DUNEPRISMTOOLSROOT}/fcl/${CONFIG_FCL_ND}"
    exit 8
  fi
  fhicl-dump ${DUNEPRISMTOOLSROOT}/fcl/${CONFIG_FCL_ND} > ND_build_flux.fcl
fi

if [ ${DO_FD} == "1" ]; then
  if [ ! -e ${DUNEPRISMTOOLSROOT}/fcl/${CONFIG_FCL_FD} ]; then
    echo "[ERROR]: Expected to find ${DUNEPRISMTOOLSROOT}/fcl/${CONFIG_FCL_FD}"
    exit 8
  fi
  fhicl-dump ${DUNEPRISMTOOLSROOT}/fcl/${CONFIG_FCL_FD} > FD_build_flux.fcl
fi

tar -zcvf apps.${DUNEPRISMTOOLS_VERSION}.tar.gz dp_BuildFluxes libTH2Jagged.so inputs.list *build_flux.fcl

echo "jobsub_submit --group=${EXPERIMENT} --jobid-output-only --resource-provides=usage_model=OPPORTUNISTIC,OFFSITE --expected-lifetime=${LIFETIME_EXP} --disk=${DISK_EXP} -N ${NJOBSTORUN} --memory=${MEM_EXP} --cpu=1 --OS=SL6 --tar_file_name=dropbox://apps.${DUNEPRISMTOOLS_VERSION}.tar.gz file://${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/BuildFluxJob.sh ${NPPFXU} ${USE_PPFX_COMPONENTS} --only-pdg ${ONLY_PDG} ${NPERJOB} ${ND_PATH_APPEND} ${FD_PATH_APPEND}"
JID=$(jobsub_submit --group=${EXPERIMENT} --jobid-output-only --resource-provides=usage_model=OPPORTUNISTIC,OFFSITE --expected-lifetime=${LIFETIME_EXP} --disk=${DISK_EXP} -N ${NJOBSTORUN} --memory=${MEM_EXP} --cpu=1 --OS=SL6 --tar_file_name=dropbox://apps.${DUNEPRISMTOOLS_VERSION}.tar.gz file://${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/BuildFluxJob.sh ${NPPFXU} ${USE_PPFX_COMPONENTS} --only-pdg ${ONLY_PDG} ${NPERJOB} ${ND_PATH_APPEND} ${FD_PATH_APPEND})

JID=$(echo ${JID} | tr -d "\n" | tr -d " " | tr -d "\t")

cd ../
rm -rf sub_tmp

echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; jobsub_q --user=${USER} --group=${EXPERIMENT} --jobid=${JID}" > JID_${JID}.q.sh
echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; jobsub_rm --user=${USER} --group=${EXPERIMENT} --jobid=${JID}" > JID_${JID}.rm.sh
echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; mkdir ${JID}_log; cd ${JID}_log; jobsub_fetchlog --user=${USER} --group=${EXPERIMENT} --jobid=${JID}; tar -zxvf *.tgz" > JID_${JID}.fetchlog.sh

chmod +x JID_${JID}.q.sh JID_${JID}.rm.sh JID_${JID}.fetchlog.sh
