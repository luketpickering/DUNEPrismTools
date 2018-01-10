#!/bin/sh

DET_DIST_CM=57400
# FD_cm 128700000
PNFS_PATH_APPEND=""
DK2NU_INPUT_DIR=""
RUN_PLAN_XML="${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml"
NMAXJOBS=""
NPERJOB="10"
INPUTLIST=""
LIFETIME_EXP="30m"
DISK_EXP="1GB"
MEM_EXP="2GB"
BINNING_DESCRIPTOR="0,0.5,1_3:0.25,3_4:0.5,4_10:1,10_20:2"
REUSEPARENTS="1"
SPECARG="0"

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

      -d|--detector-distance)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      DET_DIST_CM="$2"
      echo "[OPT]: Detector assumed to be at z=${DET_DIST_CM} cm in beam simulation coordinates."
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

      -r|--runplan)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      RUN_PLAN_XML=$(readlink -f $2)
      echo "[OPT]: Using runplan xml file: \"${RUN_PLAN_XML}\"."
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

      -N|--N-max-jobs)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      NMAXJOBS="$2"
      echo "[OPT]: Running a maximum of: \"${NMAXJOBS}\"."
      shift # past argument
      ;;

      -b|--binning)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      BINNING_DESCRIPTOR="$2"
      echo "[OPT]: Using binning descriptor: \"${BINNING_DESCRIPTOR}\"."
      shift # past argument
      ;;

      -S|--species)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      SPECARG="$2"
      echo "[OPT]: Only running for species, PDG = \"${SPECARG}\"."
      shift # past argument
      ;;

      -P|--no-reuse-parents)

      REUSEPARENTS="0"
      echo "[OPT]: Will use each decay parent once."
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
      echo -e "\t-r|--runplan               : Run Plan XML file describing detector stops."
      echo -e "\t-n|--n-per-job             : Number of files to run per job. (default: 10)"
      echo -e "\t-b|--binning               : dp_BuildFluxes variable binning descriptor. (default: 0,0.5,1_3:0.25,3_4:0.5,4_10:1,10_20:2)"
      echo -e "\t-N|--NMAXJobs              : Maximum number of jobs to submit."
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

NJOBSTORUN=$(python -c "from math import ceil; print int(ceil(float(${NINPUTS})/float(${NPERJOB}))); ")
if [ ! -z ${NMAXJOBS} ]; then
  NJOBSTORUN=$(python -c "print min(${NMAXJOBS},${NJOBSTORUN})")
fi

echo "[INFO]: Found ${NINPUTS} inputs, will run ${NJOBSTORUN} jobs."

if [ ! "${RUN_PLAN_XML}" ] || [ ! -e "${RUN_PLAN_XML}" ]; then
  echo "[ERROR]: No or non-existant runplan file: \"${RUN_PLAN_XML}\"."
  exit 1
fi

if [ ! -e ${DUNEPRISMTOOLSROOT}/bin/dp_BuildFluxes ]; then
  echo "[ERROR]: It appears that DUNEPrismTools was not built, expected to find \"${DUNEPRISMTOOLSROOT}/bin/dp_BuildFluxes\"."
  exit 1
fi

ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/flux

if [ $? -ne 0 ]; then
  mkdir -p /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/flux
  ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/flux
  if [ $? -ne 0 ]; then
    echo "Unable to make /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/flux."
    exit 7
  fi
fi

ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/logs

if [ $? -ne 0 ]; then
  mkdir -p /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/logs
  ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/logs
  if [ $? -ne 0 ]; then
    echo "Unable to make /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/logs."
    exit 7
  fi
fi

cp ${DUNEPRISMTOOLSROOT}/bin/dp_* .
cp ${RUN_PLAN_XML} ./runplan.xml

tar -zcvf apps.@DUNEPrismTools_VERSION_STRING@.tar.gz dp_* inputs.list runplan.xml

JID=$(jobsub_submit --group=${EXPERIMENT} --jobid-output-only --resource-provides=usage_model=OPPORTUNISTIC --expected-lifetime=${LIFETIME_EXP} --disk=${DISK_EXP} -N ${NJOBSTORUN} --memory=${MEM_EXP} --cpu=1 --OS=SL6 --tar_file_name=dropbox://apps.@DUNEPrismTools_VERSION_STRING@.tar.gz file://${DUNEPRISMTOOLSROOT}/scripts/BuildFluxJob.sh ${DET_DIST_CM} ${BINNING_DESCRIPTOR} ${REUSEPARENTS} ${SPECARG} ${PNFS_PATH_APPEND} ${NPERJOB})

cd ../
rm -r sub_tmp

echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; jobsub_q --jobid=${JID}" > JID_${JID}.q.sh
echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; jobsub_rm --jobid=${JID}" > JID_${JID}.rm.sh
echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; mkdir ${JID}_log; cd ${JID}_log; jobsub_fetchlog --jobid=${JID}; tar -zxvf *.tgz" > JID_${JID}.fetchlog.sh

chmod +x JID_${JID}.q.sh JID_${JID}.rm.sh JID_${JID}.fetchlog.sh
