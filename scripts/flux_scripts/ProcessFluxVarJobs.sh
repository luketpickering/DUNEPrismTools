# !/bin/bash

# #alignment
# for i in nu nubar; do
#   for j in DecayPipe Horn1 Horn2; do
#     for k in X Y; do
for i in nu; do
  for j in Horn1; do
    for k in X; do
      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
         --expected-walltime 2h --expected-disk 2GB \
         --expected-mem 512MB -F build_ND_fluxes_${i}mode.fcl \
         -p Alignment/DUNEPrismFluxes/ND_${i}/${j}${k}/uncert_binning \
         -i /pnfs/dune/persistent/users/picker24/Alignment/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${j}${k}/${i}/dk2nulite \
         -n 20 -f
    done
  done
done


# #nu
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#    --expected-walltime 2h --expected-disk 2GB \
#    --expected-mem 512MB -F build_ND_fluxes_numode_ppfx.fcl \
#    -p nominal_5E8POT_wppfx/DUNEPrismFluxes/ND_nu/uncert_binning \
#    -i /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/nu/dk2nulite \
#    -n 20 -f
#
# #nubar
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#    --expected-walltime 2h --expected-disk 2GB \
#    --expected-mem 512MB -F build_ND_fluxes_nubarmode_ppfx.fcl \
#    -p nominal_5E8POT_wppfx/DUNEPrismFluxes/ND_nubar/uncert_binning  \
#    -i /pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/nubar/dk2nulite \
#    -n 20 -f

#With focussing
# for i in nu nubar; do
#   for j in p1; do
#     for k in WL HC DPR; do
#       ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#          --expected-walltime 2h --expected-disk 2GB \
#          --expected-mem 512MB -F build_ND_fluxes_${i}mode.fcl \
#          -p Focussing/DUNEPrismFluxes/ND_${i}/${k}${j}/uncert_binning \
#          -i /pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite \
#          -n 20 -f
#     done
#   done
# done
