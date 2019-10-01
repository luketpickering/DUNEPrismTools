#!/bin/sh

PNFS_PATH_APPEND=""
NJOBSTORUN="0"
LIFETIME_EXP="3h"
DISK_EXP="1GB"
MEM_EXP="1GB"
USE_PPFX="0"
USE_ALL_PPFX_WEIGHTS="0"
MACRO=""
NPROTONS="100000"
JOBNAME="RunG4LBNE_PPFX_Makedk2nuLite"

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

      -P|--PPFX)

      USE_PPFX="1"
      echo "[OPT]: Attempting to run PPFX."
      ;;

      --Use-All-PPFX-Weights)
      USE_ALL_PPFX_WEIGHTS="1"
      echo "[OPT]: Attempting to keep all PPFX weight branches."
      ;;

      -m|--macro)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      MACRO="$2"
      echo "[OPT]: Running the beam simulation with ${MACRO}"
      shift # past argument
      ;;

      -N|--number-of-jobs)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      NJOBSTORUN="$2"
      echo "[OPT]: Will submit \"${NJOBSTORUN}\" jobs."
      shift # past argument
      ;;

      -n|--number-of-protons)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      NPROTONS="$2"
      echo "[OPT]: Will simulate \"${NPROTONS}\" events per job."
      shift # past argument
      ;;

      --jobname)

      if [[ ${#} -lt 2 ]]; then
        echo "[ERROR]: ${1} expected a value."
        exit 1
      fi

      JOBNAME="$2"
      echo "[OPT]: Naming job \"${JOBNAME}\"."
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

      -D|--dry-run)
        echo -e "[OPT]: Performing dry run."
        DRY_RUN="1"
      ;;

      -?|--help)

      echo "[RUNLIKE] ${SCRIPTNAME}"
      echo -e "\t-m|--macro                 : g4lbnf macro to use in generation."
      echo -e "\t-p|--pnfs-path-append      : Path to append to output path: /pnfs/dune/persistent/users/${USER}/"
      echo -e "\t-P|--PPFX                  : Will run PPFX."
      echo -e "\t--Use-All-PPFX-Weights     : Keep all hadron production weight branches separate, will increase disk space usage of output file by 5-10x."
      echo -e "\t-N|--number-of-jobs        : Number of flux jobs to run."
      echo -e "\t-n|--number-of-protons     : Number of POT to throw (Default: 100000)"
      echo -e "\t--jobname                  : Rename the submitted script to make job identification simpler."
      echo -e "\t--expected-disk            : Expected disk usage to pass to jobsub. "
      echo -e "\t--expected-mem             : Expected mem usage to pass to jobsub."
      echo -e "\t--expected-walltime        : Expected disk usage to pass to jobsub."
      echo -e "\t-?|--help                  : Print this message."
      exit 0
      ;;

      *)
              # unknown option
      echo "[ERROR]: Unknown option $1"
      exit 1
      ;;
  esac
  shift # past argument or value
done

if [ "${NJOBSTORUN}" == "0" ]; then
  echo "[ERROR]: Must specify the number of jobs to run."
  exit 1
fi

if [ -e sub_tmp ]; then rm -r sub_tmp; fi

mkdir sub_tmp
cd sub_tmp

if [ -z ${MACRO} ] || [ ! -e ${MACRO} ]; then
  echo "[ERROR]: Was not passed a valid -m option: \"${MACRO}\"."
  exit 1
fi
MACRO_FILE_NAME=${MACRO##*/}

if [ ! "${EXPERIMENT}" ]; then
  echo "[ERROR]: No experiment. Expected \${EXPERIMENT} to be defined."
  exit 1
fi

if [ "${EXPERIMENT}" != "dune" ]; then
  echo "[WARN]: \${EXPERIMENT} != \"dune\", \${EXPERIMENT} = \"${EXPERIMENT}\"."
fi

if [ -z ${DUNEPRISMTOOLSROOT} ]; then
  echo "[ERROR]: \$DUNEPRISMTOOLSROOT was not defined, you must set up DUNEPrismTools before deploying."
  exit 1
fi

if [ -z ${G4LBNE_DIR} ]; then
  echo "[ERROR]: \$G4LBNE_DIR was not defined, you must set up g4lbnf before deploying."
  exit 1
fi

if [ "${USE_PPFX}" == "1" ]; then
  if [ -z ${PPFX_DIR} ]; then
    echo "[ERROR]: \$PPFX_DIR was not defined, you must set up PPFX before deploying."
    exit 1
  fi
fi

mkdir tar_scratch
cd tar_scratch

cp ${MACRO} ${MACRO_FILE_NAME}.inp
cat ${MACRO_FILE_NAME}.inp | sed "s/__NEVENTS__/${NPROTONS}/g" > ${MACRO_FILE_NAME}

echo "====START====MACRO============"
cat ${MACRO_FILE_NAME}
echo "=====END=====MACRO============"

mkdir setups
cp ${G4LBNE_DIR}/setups/setup_g4lbne_cvmfs.sh setups/setup.sh
cp -r ${G4LBNE_DIR}/{g4lbnf,libmyVersion.so,libg4lbnfDict.so,g4lbnfCint_rdict.pcm,locations,gdml} ./
cp ${DUNEPRISMTOOLSROOT}/bin/dp_MakeLitedk2nu ./
cp ${DUNEPRISMTOOLSROOT}/lib/libTH2Jagged.so ./

if [ "${USE_PPFX}" == "1"  ]; then
  mkdir ppfx
  cp -r ${PPFX_DIR}/{data,uncertainties} ppfx/
  cp ${PPFX_DIR}/scripts/inputs_default.xml ./ppfx.xml

  cp ${PPFX_DIR}/bin/Make_dk2nu_FriendTree ./
  cp ${PPFX_DIR}/lib/libppfx.so ./

  tar -zcvf g4blne_makedk2nulite.tar.gz \
    ${MACRO_FILE_NAME} \
    setups/setup.sh \
    g4lbnf libmyVersion.so libg4lbnfDict.so g4lbnfCint_rdict.pcm \
    locations/* gdml/* \
    dp_MakeLitedk2nu libTH2Jagged.so \
    Make_dk2nu_FriendTree \
    libppfx.so \
    ppfx/* ppfx.xml
else
  tar -zcvf g4blne_makedk2nulite.tar.gz \
    ${MACRO_FILE_NAME} \
    setups/setup.sh \
    g4lbnf libmyVersion.so libg4lbnfDict.so g4lbnfCint_rdict.pcm \
    locations/* gdml/* \
    dp_MakeLitedk2nu libTH2Jagged.so
fi

mv g4blne_makedk2nulite.tar.gz ../
cd ../

rm -r tar_scratch

source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh

setup jobsub_client
setup ifdhc

for dir in dk2nulite logs; do
  ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/${dir}

  if [ $? -ne 0 ]; then
    echo "Attempting to make /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/${dir}."
    mkdir -p /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/${dir}
    ifdh ls /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/${dir}
    if [ $? -ne 0 ]; then
      echo "Unable to make /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/${dir}."
      exit 7
    fi
    echo "Made /pnfs/dune/persistent/users/${USER}/${PNFS_PATH_APPEND}/${dir}."
  fi
done

cp ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/RunG4LBNE_PPFX_Makedk2nuLite.sh ./${JOBNAME}.sh

JOBSCRIPT=$(readlink -f ./${JOBNAME}.sh)

if [ -z ${JOBSCRIPT} ] || [ ! -e ${JOBSCRIPT} ]; then
  echo "[ERROR]: No such file \"${JOBSCRIPT}\". Cannot proceed"
  exit 1
fi

echo "Submitting job: jobsub_submit --group=${EXPERIMENT} --jobid-output-only --resource-provides=usage_model=OPPORTUNISTIC,OFFSITE --expected-lifetime=${LIFETIME_EXP} --disk=${DISK_EXP} -N ${NJOBSTORUN} --memory=${MEM_EXP} --cpu=1 --OS=SL7 --tar_file_name=dropbox://g4blne_makedk2nulite.tar.gz file://${JOBSCRIPT} ${MACRO_FILE_NAME} ${PNFS_PATH_APPEND} ${USE_PPFX} ${USE_ALL_PPFX_WEIGHTS}"
JID=$(jobsub_submit --group=${EXPERIMENT} --jobid-output-only --resource-provides=usage_model=OPPORTUNISTIC,OFFSITE --expected-lifetime=${LIFETIME_EXP} --disk=${DISK_EXP} -N ${NJOBSTORUN} --memory=${MEM_EXP} --cpu=1 --OS=SL7 --tar_file_name=dropbox://g4blne_makedk2nulite.tar.gz file://${JOBSCRIPT} ${MACRO_FILE_NAME} ${PNFS_PATH_APPEND} ${USE_PPFX} ${USE_ALL_PPFX_WEIGHTS})



JID=$(echo ${JID} | tr -d "\n" | tr -d " " | tr -d "\t" | sed "s/\.0//g")

echo "Submitted job with ID: ${JID}"

cd ../
rm -r sub_tmp

echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; jobsub_q --user=${USER} --group=${EXPERIMENT} --jobid=${JID}" > JID_${JID}.q.sh
echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; jobsub_rm --user=${USER} --group=${EXPERIMENT} --jobid=${JID}" > JID_${JID}.rm.sh
echo -e "#!/bin/sh\n source /cvmfs/fermilab.opensciencegrid.org/products/common/etc/setups.sh; setup jobsub_client; mkdir ${JID}_log; cd ${JID}_log; jobsub_fetchlog --user=${USER} --group=${EXPERIMENT} --jobid=${JID}; tar -zxvf *.tgz" > JID_${JID}.fetchlog.sh

chmod +x JID_${JID}.q.sh JID_${JID}.rm.sh JID_${JID}.fetchlog.sh

#to reduce job submission rate
echo "Sleeping for 3 minutes..."
sleep 2.5m
