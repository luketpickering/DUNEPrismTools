#!/bin/bash

BASEFCL=${DUNEPRISMTOOLSROOT}/fcl/flux_uncertainty_generation_comp.toconfig.fcl
IFILE=DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_uncert_jagged_opt.root
#
# for twk in  DecayPipeRadius \
#             WaterLayer \
#             HornCurrent \
#             TargetDensity \
#             Horn1XShift \
#             Horn2XShift \
#             Horn1YShift \
#             Horn2YShift \
#             BeamTheta \
#             BeamThetaPhi \
#             BeamSigma \
#             BeamOffsetX \
#             PPFX \
#             POTCounting; do
#
#   cat ${BASEFCL} | sed "s/__PROPNAME__/${twk}/g" | sed "s/__TWEAKS__/@local::${twk}/g" | sed "s/__SPECIES__/[numu,numubar,nue,nuebar]/g" | sed "s/__CONFIGS__/[ND_nu, ND_nubar, FD_nu, FD_nubar]/g" > ${twk}_.fcl
#   fhicl-dump ./${twk}_.fcl | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > ${twk}.fcl
#   dp_BuildUncertaintyMatrix --fhicl ./${twk}.fcl
#
# done
#
# for set in  FocussingTweaks \
#             AlignmentTweaks \
#             BeamAlignmentTweaks; do
#
#   cat ${BASEFCL} | sed "s/__PROPNAME__/${set}/g" | sed "s/__TWEAKS__/@sequence::${set}/g" | sed "s/__SPECIES__/[numu,numubar,nue,nuebar]/g" | sed "s/__CONFIGS__/[ND_nu, ND_nubar, FD_nu, FD_nubar]/g" > ${set}_.fcl
#   fhicl-dump ./${set}_.fcl | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > ${set}.fcl
#   dp_BuildUncertaintyMatrix --fhicl ./${set}.fcl
#
# done
#
# for twk in  PPFX_abs\
#             PPFX_att\
#             PPFX_ttpCpi\
#             PPFX_ttpCk\
#             PPFX_ttnCpi\
#             PPFX_ttpCnu\
#             PPFX_ttnua\
#             PPFX_ttmesinc\
#             PPFX_oth; do
#
#   cat ${BASEFCL} | sed "s/__PROPNAME__/${twk}_nu/g" | sed "s/__TWEAKS__/@local::${twk}/g" | sed "s/__SPECIES__/[numu]/g" | sed "s/__CONFIGS__/[ND_nu,FD_nu]/g" > ${twk}_nu_.fcl
#   fhicl-dump ./${twk}_nu_.fcl | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > ${twk}_nu.fcl
#   dp_BuildUncertaintyMatrix --fhicl ./${twk}_nu.fcl
#
#   cat ${BASEFCL} | sed "s/__PROPNAME__/${twk}_nubar/g" | sed "s/__TWEAKS__/@local::${twk}/g" | sed "s/__SPECIES__/[numubar]/g" | sed "s/__CONFIGS__/[ND_nubar,FD_nubar]/g" > ${twk}_nubar_.fcl
#   fhicl-dump ./${twk}_nubar_.fcl | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > ${twk}_nubar.fcl
#   dp_BuildUncertaintyMatrix --fhicl ./${twk}_nubar.fcl
#
# done
#
# cat ${BASEFCL} | sed "s/__PROPNAME__/Total/g" | sed "s/__TWEAKS__/@local::PPFX,@local::POTCounting,@sequence::FocussingTweaks,@sequence::AlignmentTweaks,@sequence::BeamAlignmentTweaks/g" | sed "s/__SPECIES__/[numu,numubar,nue,nuebar]/g" | sed "s/__CONFIGS__/[ND_nu, ND_nubar, FD_nu, FD_nubar]/g" > Total_.fcl
# fhicl-dump ./Total_.fcl | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > Total.fcl
# dp_BuildUncertaintyMatrix --fhicl ./Total.fcl
#
# cat ${BASEFCL} | sed "s/__PROPNAME__/Total_onaxis/g" | sed "s/__TWEAKS__/@local::PPFX,@local::POTCounting,@sequence::FocussingTweaks,@sequence::AlignmentTweaks,@sequence::BeamAlignmentTweaks/g" | sed "s/__SPECIES__/[numu,numubar,nue,nuebar]/g" | sed "s/__CONFIGS__/[ND_nu, ND_nubar, FD_nu, FD_nubar]/g" | sed "s/# ND_nu_OffAxisBin/ND_nu_OffAxisBin/g" | sed "s/# ND_nubar_OffAxisBin/ND_nubar_OffAxisBin/g" > Total_onaxis_.fcl
# fhicl-dump ./Total_onaxis_.fcl | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > Total_onaxis.fcl
# dp_BuildUncertaintyMatrix --fhicl ./Total_onaxis.fcl
#

# cat ${BASEFCL} | sed "s/__PROPNAME__/TotalNonHP/g" | sed "s/__TWEAKS__/@local::POTCounting,@sequence::FocussingTweaks,@sequence::AlignmentTweaks,@sequence::BeamAlignmentTweaks/g" | sed "s/__SPECIES__/[numu,numubar,nue,nuebar]/g" | sed "s/__CONFIGS__/[ND_nu, ND_nubar, FD_nu, FD_nubar]/g" > TotalNonHP_.fcl
# fhicl-dump ./TotalNonHP_.fcl | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > TotalNonHP.fcl
# dp_BuildUncertaintyMatrix --fhicl ./TotalNonHP.fcl

# rm *.fcl
#
# fhicl-dump flux_ratio_plots.fcl | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > flux_ratio_plots.fcl
#
# dp_NDFDFluxRatioPlots ./flux_ratio_plots.fcl
#
# rm flux_ratio_plots.fcl

# fhicl-dump flux_ratio_plots_ppfxall_nu.fcl | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > flux_ratio_plots_ppfxall_nu.fcl

# dp_NDFDFluxRatioPlots ./flux_ratio_plots_ppfxall_nu.fcl

# rm flux_ratio_plots_ppfxall_nu.fcl
#
# fhicl-dump flux_ratio_plots_ppfxall_nubar.fcl | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > flux_ratio_plots_ppfxall_nubar.fcl
#
# dp_NDFDFluxRatioPlots ./flux_ratio_plots_ppfxall_nubar.fcl
#
# rm flux_ratio_plots_ppfxall_nubar.fcl

# fhicl-dump flux_ratio_plots_nonHP.fcl | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > flux_ratio_plots_nonHP.fcl
#
# dp_NDFDFluxRatioPlots ./flux_ratio_plots_nonHP.fcl
#
# rm flux_ratio_plots_nonHP.fcl
