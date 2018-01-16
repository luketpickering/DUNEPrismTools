#!/bin/bash

#Uniform binning for fits

dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/nominal_1E9/DUNEPrismFluxes/FHC/uniform_binning/flux/Fluxes.*.root" \
  -o nominal.FHC.uniform.root

dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml \
  -i "/pnfs/dune/persistent/users/picker24/nominal_5E8_FD/DUNEPrismFluxes/FHC/uniform_binning/flux/Fluxes.*.root"  \
  -o nominal.FHC.FarDet.uniform.root


dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/nominal_1E9/DUNEPrismFluxes/RHC/uniform_binning/flux/Fluxes.*.root" \
  -o nominal.RHC.uniform.root

dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.LArDUNEFV_25.8mx23.8mx56.6m.xml \
  -i "/pnfs/dune/persistent/users/picker24/nominal_5E8_FD/DUNEPrismFluxes/RHC/uniform_binning/flux/Fluxes.*.root"  \
  -o nominal.RHC.FarDet.uniform.root


#Syst binning


dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/nominal_1E9/DUNEPrismFluxes/FHC/syst_binning/flux/Fluxes.*.root" \
  -o nominal.FHC.syst_binning.root

dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/nominal_1E9/DUNEPrismFluxes/RHC/syst_binning/flux/Fluxes.*.root" \
  -o nominal.RHC.syst_binning.root


# water layer
dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/FHC/syst_binning/flux/Fluxes.*.root" \
  -o WaterLayer_p1.FHC.syst_binning.root

dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/WaterLayer_p1_5E8/DUNEPrismFluxes/RHC/syst_binning/flux/Fluxes.*.root" \
  -o WaterLayer_p1.RHC.syst_binning.root

  dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/FHC/syst_binning/flux/Fluxes.*.root" \
  -o WaterLayer_m1.FHC.syst_binning.root

dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/WaterLayer_m1_5E8/DUNEPrismFluxes/RHC/syst_binning/flux/Fluxes.*.root" \
  -o WaterLayer_m1.RHC.syst_binning.root


# HC
dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/FHC/syst_binning/flux/Fluxes.*.root" \
  -o HC_p1.FHC.syst_binning.root

dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/HC_p1_5E8/DUNEPrismFluxes/RHC/syst_binning/flux/Fluxes.*.root" \
  -o HC_p1.RHC.syst_binning.root

  dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/FHC/syst_binning/flux/Fluxes.*.root" \
  -o HC_m1.FHC.syst_binning.root

dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/HC_m1_5E8/DUNEPrismFluxes/RHC/syst_binning/flux/Fluxes.*.root" \
  -o HC_m1.RHC.syst_binning.root

# DecayPipeR
dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/FHC/syst_binning/flux/Fluxes.*.root" \
  -o DecayPipeR_p1.FHC.syst_binning.root

dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_p1_5E8/DUNEPrismFluxes/RHC/syst_binning/flux/Fluxes.*.root" \
  -o DecayPipeR_p1.RHC.syst_binning.root

  dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/FHC/syst_binning/flux/Fluxes.*.root" \
  -o DecayPipeR_m1.FHC.syst_binning.root

dp_CombineBuiltFluxes -r ${DUNEPRISMTOOLSROOT}/configs/RunPlan.39mLAr.3mFV.10cm.2mx4m.xml \
  -i "/pnfs/dune/persistent/users/picker24/DecayPipeR_m1_5E8/DUNEPrismFluxes/RHC/syst_binning/flux/Fluxes.*.root" \
  -o DecayPipeR_m1.RHC.syst_binning.root
