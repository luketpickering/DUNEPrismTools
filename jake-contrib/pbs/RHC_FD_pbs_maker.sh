#!/bin/bash

FILE=~/pbs_header_4_hr

if [ ! -d "$1" ]; then
  mkdir /home/calcuttj/DUNEPrismSim/RHC_FD/$1
  echo "Made directory /home/calcuttj/DUNEPrismSim/RHC_FD/$1/"
  mkdir RHC_FD/$1
  echo "Made directory RHC_FD/$1/"
fi

cat RHC_FD_input.list | while read line 
do
  echo $line
  string=$line
  IFS='.'
  array=( $string )
  run=$(echo ${array[2]}.${array[3]})
  IFS=''
  echo $run

  OUTFILE=RHC_FD/$1/dpa.RHC_FD.${array[2]}.${array[3]}.pbs

  cat $FILE > $OUTFILE  
  echo \#PBS -l file=20gb >> ${OUTFILE}
  echo \#PBS -N ${run}_analysis >> ${OUTFILE} 
  echo DATAPATH=/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker/RHC_FD >> ${OUTFILE}
  echo XMLPATH=/mnt/research/NuInt/DUNEPrismTools/jake-contrib/${2} >> ${OUTFILE}
  echo PYEXEC=/mnt/research/NuInt/argon_box_repo/argon_box/argon_box_genie.py >> ${OUTFILE}

  echo source /mnt/research/NuInt/argon_box_repo/argon_box/setup.sh 4.10.00.p04 64 >> ${OUTFILE}
  echo export PATH=/mnt/research/NuInt/DUNEPrismTools/jake-contrib/:\$PATH >> ${OUTFILE}
  echo module load GNU/4.9 >> ${OUTFILE}

  echo cp \$\{DATAPATH\}/${line} \$TMPDIR >> ${OUTFILE}
  echo cd \$TMPDIR >> ${OUTFILE}
  

  echo PYOUTPUT=${run}.\$\{PBS_JOBID\}.argon_box.root >> ${OUTFILE}
  echo python \$PYEXEC --nevents=0 --source=$line --output=\$\{PYOUTPUT\} --enable_edepsim --detX=50. --detY=50. --detZ=200. --shift=0 >> ${OUTFILE}

  echo cp \$\{PYOUTPUT\} /mnt/scratch/calcuttj/DunePrism/RHC_FD/ >> ${OUTFILE}

  echo DunePrismAnalyzer -i \$\{PYOUTPUT\} -o ${run}.\$\{PBS_JOBID\}.DPA.root -x \$\{XMLPATH\} -n -1 >> ${OUTFILE}

  echo cp ${run}.\$\{PBS_JOBID\}.DPA.root ~/DUNEPrismSim/RHC_FD/$1/ >> ${OUTFILE}

 
done
