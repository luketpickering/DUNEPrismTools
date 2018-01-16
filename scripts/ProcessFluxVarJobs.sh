# !/bin/bash

BINNING_DESCRIPTOR_FITS="0_10:0.025"
BINNING_DESCRIPTOR_UNCERTS="0,0.5,1_3:0.25,3_4:0.5,4_10:1,10_20:2"

### Fits

# ND

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p nominal_1E9/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/nominal_1E9/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p nominal_1E9/DUNEPrismFluxes/RHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/nominal_1E9/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/slimflux/ \
  -n 3 -D

# FD

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p nominal_5E8_FD/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p nominal_5E8_FD/DUNEPrismFluxes/RHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/slimflux/ \
  -n 3 -D


### Variations

# Nominal
${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p nominal_1E9/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/nominal_1E9/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p nominal_1E9/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/nominal_1E9/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/slimflux/ \
  -n 3 -D


# Water layer
${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p WaterLayer_p1_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p WaterLayer_p1_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/antineutrino/slimflux/ \
  -n 3 -D


${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p WaterLayer_m1_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p WaterLayer_m1_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/antineutrino/slimflux/ \
  -n 3 -D


# Decay Pipe R
${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p DecayPipeR_p1_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/antineutrino/slimflux/ \
  -n 3 -D


${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p DecayPipeR_m1_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/antineutrino/slimflux/ \
  -n 3 -D


# HC
${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p HC_p1_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p HC_p1_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/antineutrino/slimflux/ \
  -n 3 -D


${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p HC_m1_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p HC_m1_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/antineutrino/slimflux/ \
  -n 3 -D


#3
${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p HC_p3_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_p3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p3/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p HC_p3_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_p3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p3/antineutrino/slimflux/ \
  -n 3 -D


${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p HC_m3_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_m3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m3/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p HC_m3_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_m3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m3/antineutrino/slimflux/ \
  -n 3 -D

#5
${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p HC_p5_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_p5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p5/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p HC_p5_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_p5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p5/antineutrino/slimflux/ \
  -n 3 -D


${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p HC_m5_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_m5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m5/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
  -p HC_m5_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_m5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m5/antineutrino/slimflux/ \
  -n 3 -D


# For syst variations

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p WaterLayer_p1_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/neutrino/slimflux/ \
  -n 3 -D
${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p WaterLayer_m1_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/neutrino/slimflux/ \
  -n 3 -D
${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/neutrino/slimflux/ \
  -n 3 -D

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p HC_p1_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/neutrino/slimflux/ \
  -n 3 -D
${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p HC_m1_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/neutrino/slimflux/ \
  -n 3 -D


${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p HC_p3_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_p3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p3/neutrino/slimflux/ \
  -n 3 -D
${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p HC_m3_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_m3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m3/neutrino/slimflux/ \
  -n 3 -D


${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p HC_p5_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/pickeVr24/HC_p5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p5/neutrino/slimflux/ \
  -n 3 -D
${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh \
  --expected-walltime 60m --expected-disk 1GB \
  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
  -p HC_m5_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
  -i /pnfs/dune/persistent/users/picker24/HC_m5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m5/neutrino/slimflux/ \
  -n 3 -D
