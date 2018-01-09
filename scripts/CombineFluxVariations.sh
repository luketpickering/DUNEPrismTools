#!/bin/bash

# For Flux Fits

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/ForPRISMFits/FHC/flux/Fluxes.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/ForPRISMFits/FHC/flux/Fluxes.*.root" -o nominal.FHC.ForFits.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/ForPRISMFits/FHC/flux/Fluxes.*.root" -o HC_m1.FHC.ForFits.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/ForPRISMFits/FHC_noreuse/flux/Fluxes.*.root" -o nominal.FHC.ForFits.noreuse.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_5E8_FD/DUNEPrismFluxes/ForPRISMFits/FHC_FD/flux/Fluxes.*.root" -o nominal.FHC.ForFits.FarDet.root

##Unif bins
dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/ForPRISMFits/FHC_unif/flux/Fluxes.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/ForPRISMFits/FHC_unif/flux/Fluxes.*.root" -o nominal.FHC.ForFits.unif.root
dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_5E8_FD/DUNEPrismFluxes/FHC_FD_unif/flux/Fluxes.*.root"  -o nominal.FHC.ForFits.FarDet.unif.root


###########################

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/FHC/flux/Fluxes.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/FHC/flux/Fluxes.*.root" -o nominal.FHC.ForUncerts.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/RHC/flux/Fluxes.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/RHC/flux/Fluxes.*.root" -o nominal.RHC.root



# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.*.root" -o DecayPipeR_m1.FHC.ForUncerts.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.*.root" -o DecayPipeR_m1.RHC.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.*.root" -o DecayPipeR_p1.FHC.ForUncerts.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.*.root" -o DecayPipeR_p1.RHC.root




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.*.root" -o HC_m1.FHC.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.*.root" -o HC_m1.RHC.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.*.root" -o HC_p1.FHC.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.*.root" -o HC_p1.RHC.root




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.*.root" -o WaterLayer_m1.FHC.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.*.root" -o WaterLayer_m1.RHC.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.*.root" -o WaterLayer_p1.FHC.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.*.root" -o WaterLayer_p1.RHC.root



#### Far Det

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_5E8_FD/DUNEPrismFluxes/FHC_FD/flux/Fluxes.*.root" -o nominal.FHC.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.*.root" -o nominal.RHC.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/FHC_FD_tshort/flux/Fluxes.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/FHC_FD_tshort/flux/Fluxes.*.root" -o nominal.FHC.FarDet.tshort.root



# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.*.root" -o DecayPipeR_m1.FHC.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.*.root" -o DecayPipeR_m1.RHC.FarDet.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.*.root" -o DecayPipeR_p1.FHC.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.*.root" -o DecayPipeR_p1.RHC.FarDet.root




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.*.root" -o HC_m1.FHC.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.*.root" -o HC_m1.RHC.FarDet.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.*.root" -o HC_p1.FHC.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.*.root" -o HC_p1.RHC.FarDet.root




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.*.root" -o WaterLayer_m1.FHC.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.*.root" -o WaterLayer_m1.RHC.FarDet.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.*.root" -o WaterLayer_p1.FHC.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.*.root" -o WaterLayer_p1.RHC.FarDet.root

