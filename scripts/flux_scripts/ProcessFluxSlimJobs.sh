#!/bin/bash

## Near detector PPFX

${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
  -p nominal_5E7POT_w_ppfx/v3r5p4_beta/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino \
  -i /pnfs/dune/persistent/users/picker24/nominal_5E7POT_w_ppfx/v3r5p4_beta/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux \
  -n 15 -P --expected-disk 4GB --expected-mem 1GB --expected-walltime 10m

### Near detector

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p nominal_7.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/nominal_7.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ \
#   -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p nominal_2.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/nominal_2.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ \
#   -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p nominal_7.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/nominal_7.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/flux/ \
#   -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p nominal_2.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/nominal_2.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/flux/ \
#   -n 10


### Far detector

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ \
#   -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/flux/ \
#   -n 10



### Variations

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p DecayPipeR_m1_5E8/DUNEPrismFluxes/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/neutrino/flux/ \
#   -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p DecayPipeR_m1_5E8/DUNEPrismFluxes/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/antineutrino/flux/ \
#   -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p DecayPipeR_p1_5E8/DUNEPrismFluxes/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/neutrino/flux/ \
#   -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p DecayPipeR_p1_5E8/DUNEPrismFluxes/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/antineutrino/flux/ \
#   -n 10



# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p HC_m1_5E8/DUNEPrismFluxes/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/neutrino/flux/ \
#   -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p HC_m1_5E8/DUNEPrismFluxes/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/antineutrino/flux/ \
#   -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p HC_p1_5E8/DUNEPrismFluxes/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/neutrino/flux/ \
#   -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p HC_p1_5E8/DUNEPrismFluxes/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/antineutrino/flux/ \
#   -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p WaterLayer_m1_5E8/DUNEPrismFluxes/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/neutrino/flux/ \
#   -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p WaterLayer_m1_5E8/DUNEPrismFluxes/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/antineutrino/flux/ \
#   -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p WaterLayer_p1_5E8/DUNEPrismFluxes/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/neutrino/flux/ \
#   -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p WaterLayer_p1_5E8/DUNEPrismFluxes/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/antineutrino/flux/ \
#   -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p HC_m3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m3/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/HC_m3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m3/neutrino/flux/ \
#   -n 10
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p HC_m3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m3/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/HC_m3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m3/antineutrino/flux/ \
#   -n 10
#
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p HC_p3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p3/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/HC_p3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p3/neutrino/flux/ \
#   -n 10
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p HC_p3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p3/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/HC_p3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p3/antineutrino/flux/ \
#   -n 10
#
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p HC_m5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m5/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/HC_m5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m5/neutrino/flux/ \
#   -n 10
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p HC_m5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m5/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/HC_m5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m5/antineutrino/flux/ \
#   -n 10
#
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p HC_p5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p5/neutrino \
#   -i /pnfs/dune/persistent/users/picker24/HC_p5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p5/neutrino/flux/ \
#   -n 10
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmSlimFluxJobs.sh \
#   -p HC_p5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p5/antineutrino \
#   -i /pnfs/dune/persistent/users/picker24/HC_p5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p5/antineutrino/flux/ \
#   -n 10
