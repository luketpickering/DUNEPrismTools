# !/bin/bash

#nu
${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
   --expected-walltime 2h --expected-disk 2GB \
   --expected-mem 512MB -F build_ND_fluxes_ppfx_oldbin.fcl \
   -p nominal_5E8POT_wppfx/DUNEPrismFluxes/ND_nu/old_uncert_binning \
   -i /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/nu/dk2nulite \
   -n 20 -f

#nubar
${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
   --expected-walltime 2h --expected-disk 2GB \
   --expected-mem 512MB -F build_ND_fluxes_ppfx_oldbin.fcl \
   -p nominal_5E8POT_wppfx/DUNEPrismFluxes/ND_nubar/old_uncert_binning  \
   -i /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/nubar/dk2nulite \
   -n 20 -f

#With focussing
for i in nu nubar; do
  for j in p1 m1; do
    for k in WL HC DPR; do

      if [ ! -e /pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite ]; then
        echo "[INFO]: No input directory, skipping."
        continue
      fi

      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
         --expected-walltime 2h --expected-disk 2GB \
         --expected-mem 512MB -F build_ND_fluxes_oldbin.fcl \
         -p Focussing/DUNEPrismFluxes/ND_${i}/${k}${j}/old_uncert_binning \
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
         --expected-walltime 2h --expected-disk 2GB \
         --expected-mem 512MB -F build_ND_fluxes_oldbin.fcl \
         -p Alignment/DUNEPrismFluxes/ND_${i}/${j}${k}/old_uncert_binning \
         -i /pnfs/dune/persistent/users/picker24/Alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i}/dk2nulite \
         -n 20 -f
    done
  done
done
