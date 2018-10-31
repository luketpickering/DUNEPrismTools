# !/bin/bash

BINNING_DESCRIPTOR_FITS="0_10:0.025"
BINNING_DESCRIPTOR_UNCERTS="0,0.5,1_3:0.25,3_4:0.5,4_10:1,10_20:2"

#### With PPFX

${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
   --expected-walltime 70m --expected-disk 3GB \
   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
   -p nominal_5E7POT_w_ppfx/DUNEPrismFluxes/FHC/uniform_binning_wider -d 57400 \
   -i /pnfs/dune/persistent/users/picker24/nominal_5E7POT_w_ppfx/v3r5p4_beta/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/slimflux \
   -n 2 -D -W " -x -0.25_45.25:0.5 -h 300 " -f -X

#### With Dan

#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#    --expected-walltime 60m --expected-disk 1GB \
#    --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#    -p nominal_1E9/DUNEPrismFluxes/RHC/uniform_binning_wider -d 57400 \
#    -i /pnfs/dune/persistent/users/picker24/nominal_1E9/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/slimflux/ \
#    -n 10 -D -W " -x -0.25_45.25:0.5 -h 3 " -f
#
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#    --expected-walltime 60m --expected-disk 1GB \
#    --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#    -p nominal_1E9/DUNEPrismFluxes/FHC/uniform_binning_wider -d 57400 \
#    -i /pnfs/dune/persistent/users/picker24/nominal_1E9/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/slimflux/ \
#    -n 10 -D -W " -x -0.25_45.25:0.5 -h 3 " -f





## FD Nominal
#${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#  --expected-walltime 30m --expected-disk 1GB \
#  --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#  -p nominal_1E9/DUNEPrismFluxes/FHC_FD/uniform_binning -d 128700000 \
#  -i /pnfs/dune/persistent/users/picker24/nominal_1E9/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/slimflux/ \
#  -n 10 -D -W "-x -12.90_12.90:25.80 -h 2260" -f


## FD uncerts

# # DecayPipeR
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 30m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC_FD/uniform_binning -d 128700000 \
#   -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/neutrino/slimflux/ \
#   -n 10 -D -W "-x -12.90_12.90:25.80 -h 2260" -f
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 30m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC_FD/uniform_binning -d 128700000 \
#   -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/neutrino/slimflux/ \
#   -n 10 -D -W "-x -12.90_12.90:25.80 -h 2260" -f
#
#
# # HC
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 30m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_p1_5E8/DUNEPrismFluxes/FHC_FD/uniform_binning -d 128700000 \
#   -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/neutrino/slimflux/ \
#   -n 10 -D -W "-x -12.90_12.90:25.80 -h 2260" -f
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 30m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_m1_5E8/DUNEPrismFluxes/FHC_FD/uniform_binning -d 128700000 \
#   -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/neutrino/slimflux/ \
#   -n 10 -D -W "-x -12.90_12.90:25.80 -h 2260" -f
#
# HC 3
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 30m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_p3_5E8/DUNEPrismFluxes/FHC_FD/uniform_binning -d 128700000 \
#   -i /pnfs/dune/persistent/users/picker24/HC_p3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p3/neutrino/slimflux/ \
#   -n 10 -D -W "-x -12.90_12.90:25.80 -h 2260" -f
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 30m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_m3_5E8/DUNEPrismFluxes/FHC_FD/uniform_binning -d 128700000 \
#   -i /pnfs/dune/persistent/users/picker24/HC_m3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m3/neutrino/slimflux/ \
#   -n 10 -D -W "-x -12.90_12.90:25.80 -h 2260" -f
#
# # HC 5
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 30m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_p5_5E8/DUNEPrismFluxes/FHC_FD/uniform_binning -d 128700000 \
#   -i /pnfs/dune/persistent/users/picker24/HC_p5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p5/neutrino/slimflux/ \
#   -n 10 -D -W "-x -12.90_12.90:25.80 -h 2260" -f
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 30m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_m5_5E8/DUNEPrismFluxes/FHC_FD/uniform_binning -d 128700000 \
#   -i /pnfs/dune/persistent/users/picker24/HC_m5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m5/neutrino/slimflux/ \
#   -n 10 -D -W "-x -12.90_12.90:25.80 -h 2260" -f

# # WaterLayer
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 30m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p WaterLayer_p1_5E8/DUNEPrismFluxes/FHC_FD/uniform_binning -d 128700000 \
#   -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/neutrino/slimflux/ \
#   -n 10 -D -W "-x -12.90_12.90:25.80 -h 2260" -f
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 30m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p WaterLayer_m1_5E8/DUNEPrismFluxes/FHC_FD/uniform_binning -d 128700000 \
#   -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/neutrino/slimflux/ \
#   -n 10 -D -W "-x -12.90_12.90:25.80 -h 2260" -f

### Fits

### Separate

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 75m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p nominal_1E9/DUNEPrismFluxes/FHC_ND_NoReuse/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/nominal_1E9/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/slimflux/ \
#   -n 3 -D -S 14 -P
#
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 75m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p nominal_5E8_FD/DUNEPrismFluxes/FHC_ND_NoReuse/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/slimflux/ \
#   -n 3 -D -S 14 -P

# ND

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p nominal_1E9/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/nominal_1E9/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p nominal_1E9/DUNEPrismFluxes/RHC/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/nominal_1E9/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/slimflux/ \
#   -n 3 -D

# FD

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p nominal_5E8_FD/DUNEPrismFluxes/FHC/uniform_binning -d 128700000 \
#   -i /pnfs/dune/persistent/users/picker24/nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/slimflux/ \
#   -n 3 -D -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p nominal_5E8_FD/DUNEPrismFluxes/RHC/uniform_binning -d 128700000 \
#   -i /pnfs/dune/persistent/users/picker24/nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/slimflux/ \
#   -n 3 -D -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml


### Variations

# Nominal
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p nominal_1E9/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/nominal_1E9/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p nominal_1E9/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/nominal_1E9/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/slimflux/ \
#   -n 3 -D


# Water layer
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p WaterLayer_p1_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p WaterLayer_p1_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/antineutrino/slimflux/ \
#   -n 3 -D


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p WaterLayer_m1_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p WaterLayer_m1_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/antineutrino/slimflux/ \
#   -n 3 -D


# Decay Pipe R
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p DecayPipeR_p1_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/antineutrino/slimflux/ \
#   -n 3 -D


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p DecayPipeR_m1_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/antineutrino/slimflux/ \
#   -n 3 -D


# HC
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_p1_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_p1_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/antineutrino/slimflux/ \
#   -n 3 -D


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_m1_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_m1_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/antineutrino/slimflux/ \
#   -n 3 -D


#3
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_p3_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_p3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p3/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_p3_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_p3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p3/antineutrino/slimflux/ \
#   -n 3 -D


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_m3_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_m3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m3/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_m3_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_m3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m3/antineutrino/slimflux/ \
#   -n 3 -D

#5
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_p5_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_p5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p5/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_p5_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_p5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p5/antineutrino/slimflux/ \
#   -n 3 -D


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_m5_5E8/DUNEPrismFluxes/FHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_m5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m5/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_m5_5E8/DUNEPrismFluxes/RHC/syst_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_m5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m5/antineutrino/slimflux/ \
#   -n 3 -D


# For syst variations

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p WaterLayer_p1_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/neutrino/slimflux/ \
#   -n 3 -D
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p WaterLayer_m1_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/neutrino/slimflux/ \
#   -n 3 -D
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/neutrino/slimflux/ \
#   -n 3 -D

# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_p1_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/neutrino/slimflux/ \
#   -n 3 -D
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_m1_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/neutrino/slimflux/ \
#   -n 3 -D


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_p3_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_p3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p3/neutrino/slimflux/ \
#   -n 3 -D
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_m3_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_m3_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m3/neutrino/slimflux/ \
#   -n 3 -D


# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_p5_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/pickeVr24/HC_p5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p5/neutrino/slimflux/ \
#   -n 3 -D
# ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmBuildFluxJobs.sh \
#   --expected-walltime 60m --expected-disk 1GB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_m5_5E8/DUNEPrismFluxes/FHC/uniform_binning -d 57400 \
#   -i /pnfs/dune/persistent/users/picker24/HC_m5_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m5/neutrino/slimflux/ \
#   -n 3 -D
