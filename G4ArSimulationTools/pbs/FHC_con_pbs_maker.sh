#!/bin/bash

FILE=~/pbs_header_4_hr

if [ ! -d "$1" ]; then
  mkdir /home/calcuttj/DUNEPrismSim/FHC/$1
  echo "Made directory /home/calcuttj/DUNEPrismSim/FHC/$1/"
  mkdir FHC/$1
  echo "Made directory FHC/$1/"
fi
 
cat FHC_box_input.list | while read line 
do
  echo $line
  string=$line
  IFS='.'
  array=( $string )
  run=$(echo ${array[0]}.${array[1]})
  IFS=''
  echo $run

  OUTFILE=FHC/$1/dpa.FHC.${array[0]}.${array[1]}.pbs

  cat $FILE > $OUTFILE  
  echo \#PBS -l file=20gb >> ${OUTFILE}
  echo \#PBS -N ${run}_condenser >> ${OUTFILE} 
  echo DATAPATH=/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker/FHC >> ${OUTFILE}
  
  echo $(ls /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker/FHC/*${array[0]}.${array[1]}.root)
  echo CODE=$(ls /mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/gevgen_fnal/rootracker/FHC/*${array[0]}.${array[1]}.root) >> ${OUTFILE}

  echo PYEXEC=/mnt/research/NuInt/argon_box_repo/argon_box/argon_box_genie.py >> ${OUTFILE}

  echo source /mnt/research/NuInt/argon_box_repo/argon_box/setup.sh 4.10.00.p04 64 >> ${OUTFILE}
  echo export PATH=/mnt/research/NuInt/DUNEPrismTools/jake-contrib/:\$PATH >> ${OUTFILE}
  echo export PATH=/mnt/research/NuInt/DUNEPrismTools/Linux/bin:\$PATH >> ${OUTFILE}
  echo module load GNU/4.9 >> ${OUTFILE}

#  echo PYOUTPUT=${run}.\*.argon_box.root >> ${OUTFILE}
  echo PYOUTPUT=/mnt/scratch/calcuttj/DunePrism/FHC/${2}/${line} >>${OUTFILE}
  echo cp \$\{PYOUTPUT\} \$TMPDIR >> ${OUTFILE}
  echo cd \$TMPDIR >> ${OUTFILE}
  

#  echo python \$PYEXEC --nevents=0 --source=$line --output=\$\{PYOUTPUT\} --enable_edepsim --detX=19.5 --detY=1.5 --detZ=2.5 --shift=1800 >> ${OUTFILE}

  echo dp_DunePrismCondenser -i \$\{PYOUTPUT\} -ir \$\{CODE\}  -dmn -3800,-150,-250 -dmx 200,150,250 -fv 50,50,50 -o ${run}.\$\{PBS_JOBID\}.DPA.root -nt 1000 >> ${OUTFILE}
 
  echo cp ${run}.\$\{PBS_JOBID\}.DPA.root ~/DUNEPrismSim/FHC/$1/ >> ${OUTFILE}
 
done
