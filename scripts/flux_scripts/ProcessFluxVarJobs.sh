# !/bin/bash

#nu
${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
   --expected-walltime 2h --expected-disk 2GB \
   --expected-mem 512MB -F build_ND_fluxes_numode_ppfx.fcl \
   -p nominal_1.5E8_wppfx/DUNEPrismFluxes/ND_nu/uncert_binning \
   -i /pnfs/dune/persistent/users/picker24/nominal_1.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/neutrino/dk2nulite \
   -n 20 -f

#nubar
${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
   --expected-walltime 2h --expected-disk 2GB \
   --expected-mem 512MB -F build_ND_fluxes_nubarmode_ppfx.fcl \
   -p nominal_2.5E8_wppfx/DUNEPrismFluxes/ND_nubar/uncert_binning  \
   -i /pnfs/dune/persistent/users/picker24/nominal_2.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/antineutrino/dk2nulite \
   -n 20 -f

#With focussing
for i in nu nubar; do
  for j in p1; do
    for k in WL HC DPR; do
      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
         --expected-walltime 2h --expected-disk 2GB \
         --expected-mem 512MB -F build_ND_fluxes_${i}mode.fcl \
         -p Focussing/DUNEPrismFluxes/ND_${i}/${k}${j}/uncert_binning \
         -i /pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite \
         -n 20 -f
    done
  done
done
