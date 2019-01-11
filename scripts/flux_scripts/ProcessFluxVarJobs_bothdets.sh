# !/bin/bash

SUBDIR="old_uncert_binning"
BINSUFFIX="_oldbin"
FORCEOVERWRITE="true"

EDISK="2GB"
ETIME="2.5h"

for i in nu nubar; do

    if [ ! -e /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i}/dk2nulite ]; then
      echo "[INFO]: No input directory for ${i} wppfx, skipping."
      continue
    fi

  ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
     --expected-walltime ${ETIME} --expected-disk ${EDISK} \
     --expected-mem 512MB -FN build_ND_fluxes_ppfx${BINSUFFIX}.fcl -FF build_FD_fluxes_ppfx${BINSUFFIX}.fcl \
     -p nominal_5E8POT_wppfx/DUNEPrismFluxes/__DET___${i}/${SUBDIR} \
     -i /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i}/dk2nulite \
     -n 20 -f
done

#With focussing
for i in nu nubar; do
  for j in p1 m1; do
    for k in WL HC DPR; do

      if [ ! -e /pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite ]; then
        echo "[INFO]: No input directory, skipping."
        continue
      fi

      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
         --expected-walltime ${ETIME} --expected-disk ${EDISK} \
         --expected-mem 512MB -FN build_ND_fluxes${BINSUFFIX}.fcl -FF build_FD_fluxes${BINSUFFIC}.fcl \
         -p Focussing/DUNEPrismFluxes/__DET___${i}/${k}${j}/${SUBDIR} \
         -i /pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite \
         -n 20 -f
    done
  done
done

#Alignment
for i in nu nubar; do
  for j in Horn1 Horn2; do
    for k in X Y XNeg; do

      if [ ! -e /pnfs/dune/persistent/users/picker24/Alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i}/dk2nulite ]; then
        echo "[INFO]: No input directory, skipping."
        continue
      fi

      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
         --expected-walltime ${ETIME} --expected-disk ${EDISK} \
         --expected-mem 512MB -FN build_ND_fluxes${BINSUFFIX}.fcl -FF build_FD_fluxes${BINSUFFIC}.fcl \
         -p Alignment/DUNEPrismFluxes/__DET___${i}/${j}${k}/${SUBDIR} \
         -i /pnfs/dune/persistent/users/picker24/Alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i}/dk2nulite \
         -n 20 -f
    done
  done
done
