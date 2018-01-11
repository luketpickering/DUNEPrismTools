# !/bin/bash

BINNING_DESCRIPTOR_FITS="0_0.5:0.1,0.5_3:0.025,3_4:0.1,4_10:0.5,10_20:2"

BINNING_DESCRIPTOR_UNCERTS="0,0.5,1_3:0.25,3_4:0.5,4_10:1,10_20:2"

# rm -rf /pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/ForPRISMFits/FHC
# rm -rf /pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/ForPRISMFits/FHC
# rm -rf /pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/ForPRISMFits/FHC

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p nominal_2.5E8/DUNEPrismFluxes/ForPRISMFits/FHC -d 57400 -i /pnfs/dune/persistent/users/picker24/nominal_2.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p nominal_7.5E8/DUNEPrismFluxes/ForPRISMFits/FHC -d 57400 -i /pnfs/dune/persistent/users/picker24/nominal_7.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p HC_m1_5E8/DUNEPrismFluxes/ForPRISMFits/FHC -d 57400 -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/neutrino/flux/ -n 10


## No Re-use
# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_FITS} -P \
#   -p nominal_7.5E8/DUNEPrismFluxes/ForPRISMFits/FHC_noreuse -d 57400 -i /pnfs/dune/persistent/users/picker24/nominal_7.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_FITS} \
#   -p nominal_5E8_FD/DUNEPrismFluxes/ForPRISMFits/FHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ -n 10

##Unif
${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
  --expected-mem 1536MB -S 14 -b "0_10:0.025" \
  -p nominal_2.5E8/DUNEPrismFluxes/ForPRISMFits/FHC_unif -d 57400 -i /pnfs/dune/persistent/users/picker24/nominal_2.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ -n 10

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
  --expected-mem 1536MB -S 14 -b "0_10:0.025" \
  -p nominal_7.5E8/DUNEPrismFluxes/ForPRISMFits/FHC_unif -d 57400 -i /pnfs/dune/persistent/users/picker24/nominal_7.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ -n 10

${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
  --expected-mem 512MB -S 14 -b "0_10:0.025" \
  -p nominal_5E8_FD/DUNEPrismFluxes/FHC_FD_unif -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ -n 10


##### Near Det

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p nominal_2.5E8/DUNEPrismFluxes/FHC -d 57400 -i /pnfs/dune/persistent/users/picker24/nominal_2.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p nominal_2.5E8/DUNEPrismFluxes/RHC -d 57400 -i /pnfs/dune/persistent/users/picker24/nominal_2.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p nominal_7.5E8/DUNEPrismFluxes/FHC -d 57400 -i /pnfs/dune/persistent/users/picker24/nominal_7.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p nominal_7.5E8/DUNEPrismFluxes/RHC -d 57400 -i /pnfs/dune/persistent/users/picker24/nominal_7.5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/flux/ -n 10



# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC -d 57400 -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p DecayPipeR_m1_5E8/DUNEPrismFluxes/RHC -d 57400 -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/antineutrino/flux/ -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC -d 57400 -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p DecayPipeR_p1_5E8/DUNEPrismFluxes/RHC -d 57400 -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/antineutrino/flux/ -n 10



# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_m1_5E8/DUNEPrismFluxes/FHC -d 57400 -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_m1_5E8/DUNEPrismFluxes/RHC -d 57400 -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/antineutrino/flux/ -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_p1_5E8/DUNEPrismFluxes/FHC -d 57400 -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_p1_5E8/DUNEPrismFluxes/RHC -d 57400 -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/antineutrino/flux/ -n 10



# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p WaterLayer_m1_5E8/DUNEPrismFluxes/FHC -d 57400 -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p WaterLayer_m1_5E8/DUNEPrismFluxes/RHC -d 57400 -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/antineutrino/flux/ -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p WaterLayer_p1_5E8/DUNEPrismFluxes/FHC -d 57400 -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 1536MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p WaterLayer_p1_5E8/DUNEPrismFluxes/RHC -d 57400 -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/antineutrino/flux/ -n 10


#### Far Det

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p nominal_5E8_FD/DUNEPrismFluxes/FHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p nominal_5E8_FD/DUNEPrismFluxes/RHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/nominal_5E8_FD/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review/antineutrino/flux/ -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p DecayPipeR_m1_5E8/DUNEPrismFluxes/RHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_m1/antineutrino/flux/ -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p DecayPipeR_p1_5E8/DUNEPrismFluxes/RHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_DecayPipeR_p1/antineutrino/flux/ -n 10



# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_m1_5E8/DUNEPrismFluxes/FHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_m1_5E8/DUNEPrismFluxes/RHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/HC_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_m1/antineutrino/flux/ -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_p1_5E8/DUNEPrismFluxes/FHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p HC_p1_5E8/DUNEPrismFluxes/RHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/HC_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_HC_p1/antineutrino/flux/ -n 10



# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p WaterLayer_m1_5E8/DUNEPrismFluxes/FHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p WaterLayer_m1_5E8/DUNEPrismFluxes/RHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_m1/antineutrino/flux/ -n 10


# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p WaterLayer_p1_5E8/DUNEPrismFluxes/FHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/neutrino/flux/ -n 10

# ${DUNEPRISMTOOLSROOT}/scripts/FarmBuildFluxJobs.sh --expected-walltime 30m --expected-disk 1536MB \
#   --expected-mem 512MB -b ${BINNING_DESCRIPTOR_UNCERTS} \
#   -p WaterLayer_p1_5E8/DUNEPrismFluxes/RHC_FD -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -d 128700000 -i /pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/v3r5p4/QGSP_BERT/OptimizedEngineeredSept2017Review_WaterLayer_p1/antineutrino/flux/ -n 10
