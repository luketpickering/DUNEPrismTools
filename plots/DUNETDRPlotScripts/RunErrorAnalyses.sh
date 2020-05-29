#!/bin/bash

BASEFCL=${DUNEPRISMTOOLSROOT}/fcl/flux_uncertainty_generation_comp.toconfig.fcl
BASEFCL_ALLPCA=${DUNEPRISMTOOLSROOT}/fcl/flux_uncertainty_generation_comp.allpca.toconfig.fcl
IFILE=DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_uncertbin_offaxis_280kAOnAxis.root

# for twk in  PPFX_abs\
#             PPFX_att\
#             PPFX_ttpCpi\
#             PPFX_ttpCk\
#             PPFX_ttnCpi\
#             PPFX_ttpCnu\
#             PPFX_ttnua\
#             PPFX_ttmesinc\
#             PPFX_oth; do

#   cat ${BASEFCL} \
	# | sed "s/__PROPNAME__/${twk}_nu/g" \
	# | sed "s/__TWEAKS__/@local::${twk}/g" \
	# | sed "s/__SPECIES__/[numu]/g" \
	# | sed "s/__CONFIGS__/[ND_nu,FD_nu]/g" > ${twk}_nu_.fcl
#   fhicl-dump ./${twk}_nu_.fcl \
	# | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > ${twk}_nu.fcl
#   dp_BuildUncertaintyMatrix --fhicl ./${twk}_nu.fcl

#   cat ${BASEFCL} \
	# | sed "s/__PROPNAME__/${twk}_nubar/g" \
	# | sed "s/__TWEAKS__/@local::${twk}/g" \
	# | sed "s/__SPECIES__/[numubar]/g" \
	# | sed "s/__CONFIGS__/[ND_nubar,FD_nubar]/g" > ${twk}_nubar_.fcl
#   fhicl-dump ./${twk}_nubar_.fcl \
	# | sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > ${twk}_nubar.fcl
#   dp_BuildUncertaintyMatrix --fhicl ./${twk}_nubar.fcl

# done

CONFIGS="[ND_nu, ND_nubar, FD_nu, FD_nubar, ND_nu_280kA, ND_nubar_280kA]"

# cat ${BASEFCL} \
# 	| sed "s/__PROPNAME__/Total/g" \
# 	| sed "s/__TWEAKS__/@local::PPFX,@local::POTCounting,@sequence::FocussingTweaks,@sequence::AlignmentTweaks,@sequence::BeamAlignmentTweaks/g" \
# 	| sed "s/__SPECIES__/[numu,numubar,nue,nuebar]/g" \
# 	| sed "s/__CONFIGS__/${CONFIGS}/g" > Total_.fcl
# fhicl-dump ./Total_.fcl \
# 	| sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > Total.fcl
# dp_BuildUncertaintyMatrix --fhicl ./Total.fcl

# cat ${BASEFCL_ALLPCA} \
# 	| sed "s/__PROPNAME__/Total_allpca/g" \
# 	| sed "s/__TWEAKS__/@local::PPFX,@local::POTCounting,@sequence::FocussingTweaks,@sequence::AlignmentTweaks,@sequence::BeamAlignmentTweaks/g" \
# 	| sed "s/__SPECIES__/[numu,numubar,nue,nuebar]/g" \
# 	| sed "s/__CONFIGS__/${CONFIGS}/g" > Total_allpca_.fcl
# fhicl-dump ./Total_allpca_.fcl \
# 	| sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > Total_allpca.fcl
# dp_BuildUncertaintyMatrix --fhicl ./Total_allpca.fcl

# cat ${BASEFCL} \
# 	| sed "s/__PROPNAME__/Total_onaxis/g" \
# 	| sed "s/__TWEAKS__/@local::PPFX,@local::POTCounting,@sequence::FocussingTweaks,@sequence::AlignmentTweaks,@sequence::BeamAlignmentTweaks/g" \
# 	| sed "s/__SPECIES__/[numu,numubar,nue,nuebar]/g" \
# 	| sed "s/__CONFIGS__/${CONFIGS}/g" \
# 	| sed "s/#ONEBINGUARD//g" > Total_onaxis_.fcl
# fhicl-dump ./Total_onaxis_.fcl \
# 	| sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > Total_onaxis.fcl
# dp_BuildUncertaintyMatrix --fhicl ./Total_onaxis.fcl

cat ${BASEFCL} \
	| sed "s/__PROPNAME__/Total_onaxis_onebin/g" \
	| sed "s/__TWEAKS__/@local::PPFX,@local::POTCounting,@sequence::FocussingTweaks,@sequence::AlignmentTweaks,@sequence::BeamAlignmentTweaks/g" \
	| sed "s/__SPECIES__/[numu,numubar,nue,nuebar]/g" \
	| sed "s/__CONFIGS__/${CONFIGS}/g" \
	| sed "s/#JUSTONEBINGUARD//g" > Total_onaxis_onebin_.fcl
fhicl-dump ./Total_onaxis_onebin_.fcl \
	| sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > Total_onaxis_onebin.fcl
dp_BuildUncertaintyMatrix --fhicl ./Total_onaxis_onebin.fcl

# cat ${BASEFCL_ALLPCA} \
# 	| sed "s/__PROPNAME__/Total_onaxis_allpca/g" \
# 	| sed "s/__TWEAKS__/@local::PPFX,@local::POTCounting,@sequence::FocussingTweaks,@sequence::AlignmentTweaks,@sequence::BeamAlignmentTweaks/g" \
# 	| sed "s/__SPECIES__/[numu,numubar,nue,nuebar]/g" \
# 	| sed "s/__CONFIGS__/${CONFIGS}/g" \
# 	| sed "s/#ONEBINGUARD//g" > Total_onaxis_allpca_.fcl
# fhicl-dump ./Total_onaxis_allpca_.fcl \
# 	| sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > Total_onaxis_allpca.fcl
# dp_BuildUncertaintyMatrix --fhicl ./Total_onaxis_allpca.fcl


# cat ${BASEFCL} \
# 	| sed "s/__PROPNAME__/TotalNonHP/g" \
# 	| sed "s/__TWEAKS__/@local::POTCounting,@sequence::FocussingTweaks,@sequence::AlignmentTweaks,@sequence::BeamAlignmentTweaks/g" \
# 	| sed "s/__SPECIES__/[numu,numubar,nue,nuebar]/g" \
# 	| sed "s/__CONFIGS__/${CONFIGS}/g" > TotalNonHP_.fcl
# fhicl-dump ./TotalNonHP_.fcl \
# 	| sed "s|Fluxes/DUNE_Flux_OffAxis_Nov2017Review_syst_shifts.root|${IFILE}|g" > TotalNonHP.fcl
# dp_BuildUncertaintyMatrix --fhicl ./TotalNonHP.fcl

# dp_NDFDFluxRatioPlots flux_ratio_plots.fcl ${IFILE} FluxErrors_UncertPoints_Total.root FluxRatios.root