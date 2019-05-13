# !/bin/bash

DET="ND"
SUBDIR="old_uncert_binning_mrad"
BINSUFFIX="_oldbin"
FORCEOVERWRITE="false"

EDISK="2GB"
ETIME="2h"

for i in nu nubar; do

    if [ ! -e /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i}/dk2nulite ]; then
      echo "[INFO]: No input directory for ${i} wppfx, skipping."
      continue
    fi

    if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/DUNEPrismFluxes/${DET}_${i}/${SUBDIR} ]; then
      echo "[INFO]: Already have ${i} wppfx not reprocessing."
      continue
    fi

  ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
     --expected-walltime ${ETIME} --expected-disk ${EDISK} \
     --expected-mem 512MB -FN build_${DET}_fluxes_mrad_ppfx${BINSUFFIX}.fcl \
     -p nominal_5E8POT_wppfx/DUNEPrismFluxes/${DET}_${i}/${SUBDIR} \
     -i /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i}/dk2nulite \
     -n 20 -f

done

#With focussing
for i in nu nubar; do
  for j in p1 m1; do
    for k in WL HC DPR TargetDensity; do

      if [ ! -e /pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite ]; then
        echo "[INFO]: No input directory, skipping."
        continue
      fi

      if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/picker24/Focussing/DUNEPrismFluxes/${DET}_${i}/${k}${j}/${SUBDIR} ]; then
        echo "[INFO]: Already have ${DET}_${i}/${k}${j} not reprocessing."
        continue
      fi

      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
         --expected-walltime ${ETIME} --expected-disk ${EDISK} \
         --expected-mem 512MB -FN build_${DET}_fluxes_mrad${BINSUFFIX}.fcl \
         -p Focussing/DUNEPrismFluxes/${DET}_${i}/${k}${j}/${SUBDIR} \
         -i /pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite \
         -n 20 -f
    done
  done
done

#Alignment
for i in nu nubar; do
  for j in Horn1 Horn2; do
    for k in X Y XNeg X3mm XNeg3mm; do

      if [ ! -e /pnfs/dune/persistent/users/picker24/Alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i}/dk2nulite ]; then
        echo "[INFO]: No input directory, skipping."
        continue
      fi

      if [ ${FORCEOVERWRITE} != "true" ] && [ -e /pnfs/dune/persistent/users/picker24/Alignment/DUNEPrismFluxes/${DET}_${i}/${j}${k}/${SUBDIR} ]; then
        echo "[INFO]: Already have ${DET}_${i}/${k}${j} not reprocessing."
        continue
      fi

      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
         --expected-walltime ${ETIME} --expected-disk ${EDISK} \
         --expected-mem 512MB -FN build_${DET}_fluxes_mrad${BINSUFFIX}.fcl \
         -p Alignment/DUNEPrismFluxes/${DET}_${i}/${j}${k}/${SUBDIR} \
         -i /pnfs/dune/persistent/users/picker24/Alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i}/dk2nulite \
         -n 20 -f
    done
  done
done
