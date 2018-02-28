#!/bin/bash

FILE=~/pbs_header_1_hr

if [ ! -d "FHC/$1/processed$2/" ]; then
  mkdir FHC/$1/processed$2/
  echo "Made directory FHC/$1/processed$2"
fi
if [ ! -d "/home/calcuttj/DUNEPrismSim/FHC/$1/processed$2/" ]; then
  mkdir /home/calcuttj/DUNEPrismSim/FHC/$1/processed$2/
  echo "Made directory /home/calcuttj/DUNEPrismSim/FHC/$1/processed$2/"
fi
 
cat FHC_processor_input.list | while read line 
do
  echo $line
  string=$line
  IFS='.'
  array=( $string )
  run=$(echo ${array[0]}.${array[1]})
  IFS=''
 # echo $run

  OUTFILE=FHC/$1/processed$2/processed.FHC.${array[0]}.${array[1]}.pbs

  cat $FILE > $OUTFILE  
  echo \#PBS -l file=20gb >> ${OUTFILE}
  echo \#PBS -N ${run}_analysis >> ${OUTFILE} 
  echo DATAPATH=/home/calcuttj/DUNEPrismSim/FHC/${1}/ >> ${OUTFILE}
  echo XMLPATH=/mnt/research/NuInt/DUNEPrismTools/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml >> ${OUTFILE}

  echo source /mnt/research/NuInt/argon_box_repo/argon_box/setup.sh 4.10.00.p04 64 >> ${OUTFILE}
  echo export PATH=/mnt/research/NuInt/DUNEPrismTools/jake-contrib/:\$PATH >> ${OUTFILE}
  #echo export PATH=/mnt/research/NuInt/DUNEPrismTools/Linux/bin:\$PATH >> ${OUTFILE}
  echo module load GNU/4.9 >> ${OUTFILE}

  echo cp \$\{DATAPATH\}/${line} \$TMPDIR >> ${OUTFILE}
  echo cd \$TMPDIR >> ${OUTFILE}
  
  echo dp_FullDetTreeStopProcessor -i ${line} -o ${run}.\$\{PBS_JOBID\}.processed.root -r \$\{XMLPATH\}  >> ${OUTFILE}
 
  echo cp ${run}.\$\{PBS_JOBID\}.processed.root ~/DUNEPrismSim/FHC/$1/processed$2/ >> ${OUTFILE}

 
done    
