#!/bin/sh

PNFS_PATH_APPEND=""
DK2NU_INPUT_DIR=""
FLUX_WINDOW_NAME="DUNEPrismFluxWindow"
FLUX_WINDOW_XML="${DUNEPRISMTOOLSROOT}/configs/DUNEPrismFluxWindow.xml"
GDML_FILE="${DUNEPRISMTOOLSROOT}/configs/DUNEPrismLArBox.geom.manual.gdml"
NMAXJOBS=""
DOLOGS=""

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

      -g|--gdml)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      GDML_FILE=$(readlink -f $2)
      echo "[OPT]: Using gdml file: \"${GDML_FILE}\"."
      shift # past argument
      ;;

      -w|--flux-window-file)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      FLUX_WINDOW_XML=$(readlink -f $2)
      echo "[OPT]: Using flux window from file: \"${FLUX_WINDOW_XML}\"."
      shift # past argument
      ;;

      -f|--flux-window-name)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      FLUX_WINDOW_NAME=$(readlink -f $2)
      echo "[OPT]: Using flux window named: \"${FLUX_WINDOW_NAME}\"."
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

      -L|--log)

      DOLOGS="TRUE"
      echo "[OPT]: Writing out log files."
      ;;

      -?|--help)

      echo "[RUNLIKE] ${SCRIPTNAME}"
      echo -e "\t-p|--pnfs-path-append      : Path to append to output path: /pnfs/dune/persistent/users/${USER}/"
      echo -e "\t-i|--dk2nu-input-directory : Input directory to search for dk2nu.root files"
      echo -e "\t-g|--gdml                  : GDML geometry file"
      echo -e "\t-w|--flux-window-file      : dk2nu flux window file"
      echo -e "\t-f|--flux-window-name      : dk2nu flux window paramset name"
      echo -e "\t-N|--NMAXJobs              : Maximum number of jobs to submit."
      echo -e "\t-L|--log                   : Copy back gevgen_fnal log files."
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

if [ ! "${DK2NU_INPUT_DIR}" ] || [ ! -e "${DK2NU_INPUT_DIR}" ]; then
  echo "[ERROR]: No or non-existant dk2nu input dir: \"${DK2NU_INPUT_DIR}\"."
  exit 1
fi

PNFS_PREFIX=${DK2NU_INPUT_DIR%%/dune*}

NINPUTS=0
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

NJOBSTORUN=${NINPUTS}

if [ ! -z ${NMAXJOBS} ]; then
  NJOBSTORUN=$(python -c "print min(${NMAXJOBS},${NJOBSTORUN})")
fi

echo "[INFO]: Found ${NINPUTS} in ${DK2NU_INPUT_DIR}."
echo "[INFO]: Submitting ${NJOBSTORUN} jobs."

if [ ! "${GDML_FILE}" ] || [ ! -e "${GDML_FILE}" ]; then
  echo "[ERROR]: No or non-existant geometry file: \"${GDML_FILE}\"."
  exit 1
fi

if [ ! "${FLUX_WINDOW_XML}" ] || [ ! -e "${FLUX_WINDOW_XML}" ]; then
  echo "[ERROR]: No or non-existant flux window file: \"${FLUX_WINDOW_XML}\"."
  exit 1
fi

if [ ! "${FLUX_WINDOW_NAME}" ]; then
  echo "[ERROR]: No flux window name supplied, please specify a -f parameter."
  exit 1
fi

ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/gevgen_fnal/gntp

if [ $? -ne 0 ]; then
  mkdir -p /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/gevgen_fnal/gntp
  ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/gevgen_fnal/gntp
  if [ $? -ne 0 ]; then
    echo "Unable to make /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/gevgen_fnal/gntp."
    exit 7
  fi
fi

ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/gevgen_fnal/rootracker

if [ $? -ne 0 ]; then
  mkdir -p /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/gevgen_fnal/rootracker
  ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/gevgen_fnal/rootracker
  if [ $? -ne 0 ]; then
    echo "Unable to make /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/gevgen_fnal/rootracker."
    exit 7
  fi
fi

if [ ! -z ${DOLOGS} ]; then
  ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/gevgen_fnal/logs

  if [ $? -ne 0 ]; then
    mkdir -p /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/gevgen_fnal/logs
    ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/gevgen_fnal/logs
    if [ $? -ne 0 ]; then
      echo "Unable to make /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/gevgen_fnal/logs."
      exit 7
    fi
  fi
fi

cp ${FLUX_WINDOW_XML} ./dk2nu_FluxWindow.xml
cp ${GDML_FILE} ./geom.gdml

tar -zcvf DUNEPrism_gevgen_fnal.params.tar.gz geom.gdml dk2nu_FluxWindow.xml inputs.list

JID=$(jobsub_submit --group=${EXPERIMENT} --jobid-output-only --resource-provides=usage_model=OPPORTUNISTIC --expected-lifetime=8h --disk=1GB -N ${NJOBSTORUN} --memory=2GB --cpu=1 --OS=SL6 --tar_file_name=dropbox://DUNEPrism_gevgen_fnal.params.tar.gz file://${DUNEPRISMTOOLSROOT}/scripts/ThrowGENIEEventsJob.sh ${FLUX_WINDOW_NAME} ${PNFS_PATH_APPEND} ${DOLOGS})

cd ../
rm -r sub_tmp

echo "Submitted job: JID${JID}"

echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; jobsub_q --jobid=${JID}" > JID_${JID}.q.sh
echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; jobsub_rm --jobid=${JID}" > JID_${JID}.rm.sh
echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; mkdir ${JID}_log; cd ${JID}_log; jobsub_fetchlog --jobid=${JID}; tar -zxvf *.tgz" > JID_${JID}.fetchlog.sh

chmod +x JID_${JID}.q.sh JID_${JID}.rm.sh JID_${JID}.fetchlog.sh
