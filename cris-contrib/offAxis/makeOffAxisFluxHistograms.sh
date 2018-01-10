
#for angle_in_mrad in 0.0 1.0 2.0 3.0 4.0 5.0 10.0 20.0 30.0

MIN_MRAD=0.
MAX_MRAD=113.
STEP_MRAD=0.872665
for angle_in_mrad in 0.0
#for angle_in_mrad in `seq ${MIN_MRAD} ${STEP_MRAD} ${MAX_MRAD}` # 0.05 degree steps up to 6.45 degrees
do

# Use small angle approximation
#NEAR_DET_Z_IN_CM=57400.0;
NEAR_DET_Z_IN_CM=130000000.0; #FD
NEAR_DET_X_IN_CM=0.0
NEAR_DET_Y_IN_CM=$(echo "$NEAR_DET_Z_IN_CM*${angle_in_mrad}*.001" | bc)

echo "Calculating flux for position ($NEAR_DET_X_IN_CM,$NEAR_DET_Y_IN_CM,$NEAR_DET_Z_IN_CM)"

#root -q -b makeFluxHistograms.C\(\"$PWD\",\"100000\",\"$NEAR_DET_X_IN_CM\",\"$NEAR_DET_Y_IN_CM\",\"$NEAR_DET_Z_IN_CM\"\,\""OffAxis_${angle_in_mrad}_mrad"\"\)

# CD1 CDR reference beam w/ 120 GeV protons and 200 kA horn currents

# Split 1000 file in equal chunks for each step
N_FILES=1000
#if [ $angle_in_mrad -eq 0 ]
if (( $(echo "$angle_in_mrad == 0." |bc -l) ))
then
    I_STEP=0
else
    I_STEP=$( echo "($angle_in_mrad - $MIN_MRAD)/$STEP_MRAD" | bc -l )
fi
N_STEPS=$( echo "($MAX_MRAD - $MIN_MRAD)/$STEP_MRAD" | bc -l )
FILES_PER_STEP=$( echo "$N_FILES / $N_STEPS" | bc -l )

FIRST_FILE=`printf %04.0f $( echo "$FILES_PER_STEP * $I_STEP + 1" | bc -l )`
LAST_FILE=`printf %04.0f $( echo "$FILES_PER_STEP * ( $I_STEP + 1)" | bc -l)`

echo "I_STEP" $I_STEP "N_STEPS" $N_STEPS "FILES_PER_STEP" $FILES_PER_STEP "FIRST_FILE" $FIRST_FILE "LAST_FILE" $LAST_FILE 

echo "#!/bin/bash" > script_oa_${angle_in_mrad}.sh
echo "cd /mnt/home/f0003917/DunePRISM/offAxis"  >> script_oa_${angle_in_mrad}.sh
echo "mkdir temp_"${angle_in_mrad}  >> script_oa_${angle_in_mrad}.sh
echo "cd temp_"${angle_in_mrad}  >> script_oa_${angle_in_mrad}.sh
echo "ln -s ../*.C ." >> script_oa_${angle_in_mrad}.sh
echo "ln -s ../*.h ." >> script_oa_${angle_in_mrad}.sh
echo "ln -s ../OscLib ." >> script_oa_${angle_in_mrad}.sh
echo "ln -s ../data ." >> script_oa_${angle_in_mrad}.sh
#echo "ln -s ../g4lbne*.root ." >> script_oa_${angle_in_mrad}.sh
echo "ln -s /home/calcuttj/beam_inputs/nu/g4lbne_v3r4p2_QGSP_BERT_CP_run15_12388_80GeV_neutrino_0{$FIRST_FILE..$LAST_FILE}.dk2nu.root ." >> script_oa_${angle_in_mrad}.sh
echo "echo 'Starting flux histogram script'" >> script_oa_${angle_in_mrad}.sh
echo date >> script_oa_${angle_in_mrad}.sh
echo root -q -b makeFluxHistograms.C\'\(\"${PWD}"/temp_"${angle_in_mrad}\",\"100000\",\"$NEAR_DET_X_IN_CM\",\"$NEAR_DET_Y_IN_CM\",\"$NEAR_DET_Z_IN_CM\"\,\""OffAxis_${angle_in_mrad}_mrad"\"\)\' >> script_oa_${angle_in_mrad}.sh
echo "mv histos_*.root ../" >> script_oa_${angle_in_mrad}.sh
echo "cd - " >> script_oa_${angle_in_mrad}.sh
echo "rm -rf temp_"${angle_in_mrad} >> script_oa_${angle_in_mrad}.sh
echo "echo 'Flux histogram script ends'" >> script_oa_${angle_in_mrad}.sh
echo date >> script_oa_${angle_in_mrad}.sh


done
