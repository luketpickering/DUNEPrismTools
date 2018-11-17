# !/bin/bash

BINNING_DESCRIPTOR_FITS="0_10:0.025"
BINNING_DESCRIPTOR_UNCERTS="0,0.5,1_3:0.25,3_4:0.5,4_10:1,10_20:2"

#### With PPFX

#nu
${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
   --expected-walltime 2h --expected-disk 2GB \
   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
   -p nominal_1.5E8_wppfx/DUNEPrismFluxes/ND_nu/uncert_binning -d 57400 \
   -i /pnfs/dune/persistent/users/picker24/nominal_1.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/neutrino/dk2nulite \
   -n 20 -D -W " -x -0.25_45.25:0.5 -h 300 " -f -X

${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
   --expected-walltime 2h --expected-disk 2GB \
   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_FITS} \
   -p nominal_1.5E8_wppfx/DUNEPrismFluxes/ND_nu/fit_binning -d 57400 \
   -i /pnfs/dune/persistent/users/picker24/nominal_1.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/neutrino/dk2nulite \
   -n 20 -D -W " -x -0.25_45.25:0.5 -h 300 " -f -X

#nubar
${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
   --expected-walltime 2h --expected-disk 2GB \
   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
   -p nominal_2.5E8_wppfx/DUNEPrismFluxes/ND_nubar/uncert_binning -d 57400 \
   -i /pnfs/dune/persistent/users/picker24/nominal_2.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/antineutrino/dk2nulite \
   -n 20 -D -W " -x -0.25_45.25:0.5 -h 300 " -f -X

${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
   --expected-walltime 2h --expected-disk 2GB \
   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_FITS} \
   -p nominal_2.5E8_wppfx/DUNEPrismFluxes/ND_nubar/fit_binning -d 57400 \
   -i /pnfs/dune/persistent/users/picker24/nominal_2.5E8_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/antineutrino/dk2nulite \
   -n 20 -D -W " -x -0.25_45.25:0.5 -h 300 " -f -X

#With focussing
for i in nu nubar; do
  for j in p1; do
    for k in WL HC DPR; do
      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
         --expected-walltime 2h --expected-disk 2GB \
         --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
         -p Focussing/DUNEPrismFluxes/ND_${i}/${k}${j}/uncert_binning -d 57400 \
         -i /pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite \
         -n 20 -D -W " -x -0.25_45.25:0.5 -h 300 " -f

      ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
         --expected-walltime 2h --expected-disk 2GB \
         --expected-mem 512MB -b ${BINNING_DESCRIPTOR_FITS} \
         -p Focussing/DUNEPrismFluxes/ND_${i}/${k}${j}/fit_binning -d 57400 \
         -i /pnfs/dune/persistent/users/picker24/Focussing/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${k}${j}/${i}/dk2nulite \
         -n 20 -D -W " -x -0.25_45.25:0.5 -h 300 " -f
    done
  done
done
