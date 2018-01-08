#!/bin/bash

# For Flux Fits

dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/ForPRISMFits/FHC/flux/Fluxes.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/ForPRISMFits/FHC/flux/Fluxes.*.root" -o nominal.FHC.ForFits.root

dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/ForPRISMFits/FHC/flux/Fluxes.*.root" -o HC_m1.FHC.ForFits.root

###########################

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Coarse.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Coarse.*.root" -o nominal.FHC.Coarse.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Coarse.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Coarse.*.root" -o nominal.RHC.Coarse.root



# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Coarse.*.root" -o DecayPipeR_m1.FHC.Coarse.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Coarse.*.root" -o DecayPipeR_m1.RHC.Coarse.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Coarse.*.root" -o DecayPipeR_p1.FHC.Coarse.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Coarse.*.root" -o DecayPipeR_p1.RHC.Coarse.root




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Coarse.*.root" -o HC_m1.FHC.Coarse.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Coarse.*.root" -o HC_m1.RHC.Coarse.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Coarse.*.root" -o HC_p1.FHC.Coarse.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Coarse.*.root" -o HC_p1.RHC.Coarse.root




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Coarse.*.root" -o WaterLayer_m1.FHC.Coarse.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Coarse.*.root" -o WaterLayer_m1.RHC.Coarse.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Coarse.*.root" -o WaterLayer_p1.FHC.Coarse.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Coarse.*.root" -o WaterLayer_p1.RHC.Coarse.root



################Fine###########################



# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Fine.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Fine.*.root" -o nominal.FHC.Fine.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Fine.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Fine.*.root" -o nominal.RHC.Fine.root



# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Fine.*.root" -o DecayPipeR_m1.FHC.Fine.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Fine.*.root" -o DecayPipeR_m1.RHC.Fine.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Fine.*.root" -o DecayPipeR_p1.FHC.Fine.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Fine.*.root" -o DecayPipeR_p1.RHC.Fine.root




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Fine.*.root" -o HC_m1.FHC.Fine.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Fine.*.root" -o HC_m1.RHC.Fine.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Fine.*.root" -o HC_p1.FHC.Fine.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Fine.*.root" -o HC_p1.RHC.Fine.root




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Fine.*.root" -o WaterLayer_m1.FHC.Fine.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Fine.*.root" -o WaterLayer_m1.RHC.Fine.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Fine.*.root" -o WaterLayer_p1.FHC.Fine.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Fine.*.root" -o WaterLayer_p1.RHC.Fine.root






################Optimized###########################




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Optimized.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Optimized.*.root" -o nominal.FHC.Optimized.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Optimized.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Optimized.*.root" -o nominal.RHC.Optimized.root



# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Optimized.*.root" -o DecayPipeR_m1.FHC.Optimized.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Optimized.*.root" -o DecayPipeR_m1.RHC.Optimized.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Optimized.*.root" -o DecayPipeR_p1.FHC.Optimized.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Optimized.*.root" -o DecayPipeR_p1.RHC.Optimized.root




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Optimized.*.root" -o HC_m1.FHC.Optimized.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Optimized.*.root" -o HC_m1.RHC.Optimized.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Optimized.*.root" -o HC_p1.FHC.Optimized.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Optimized.*.root" -o HC_p1.RHC.Optimized.root




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Optimized.*.root" -o WaterLayer_m1.FHC.Optimized.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Optimized.*.root" -o WaterLayer_m1.RHC.Optimized.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/FHC/flux/Fluxes.Optimized.*.root" -o WaterLayer_p1.FHC.Optimized.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2.5x4m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/RHC/flux/Fluxes.Optimized.*.root" -o WaterLayer_p1.RHC.Optimized.root




#### Far Det

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.Optimized.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.Optimized.*.root" -o nominal.FHC.Optimized.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.Optimized.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.Optimized.*.root" -o nominal.RHC.Optimized.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/nominal_2.5E8/DUNEPrismFluxes/FHC_FD_tshort/flux/Fluxes.Optimized.*.root" -i "/pnfs/dune/persistent/users/picker24/nominal_7.5E8/DUNEPrismFluxes/FHC_FD_tshort/flux/Fluxes.Optimized.*.root" -o nominal.FHC.Optimized.FarDet.tshort.root



# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.Optimized.*.root" -o DecayPipeR_m1.FHC.Optimized.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.Optimized.*.root" -o DecayPipeR_m1.RHC.Optimized.FarDet.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.Optimized.*.root" -o DecayPipeR_p1.FHC.Optimized.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.Optimized.*.root" -o DecayPipeR_p1.RHC.Optimized.FarDet.root




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.Optimized.*.root" -o HC_m1.FHC.Optimized.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.Optimized.*.root" -o HC_m1.RHC.Optimized.FarDet.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.Optimized.*.root" -o HC_p1.FHC.Optimized.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.Optimized.*.root" -o HC_p1.RHC.Optimized.FarDet.root




# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.Optimized.*.root" -o WaterLayer_m1.FHC.Optimized.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.Optimized.*.root" -o WaterLayer_m1.RHC.Optimized.FarDet.root


# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/FHC_FD/flux/Fluxes.Optimized.*.root" -o WaterLayer_p1.FHC.Optimized.FarDet.root

# dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/RHC_FD/flux/Fluxes.Optimized.*.root" -o WaterLayer_p1.RHC.Optimized.FarDet.root

