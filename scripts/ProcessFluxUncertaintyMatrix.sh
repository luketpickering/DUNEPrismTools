#!/bin/bash

#-n Nominal.root -v Name,Up.root,Down.root -E 3 -X 35 -M 1 -o Output.root

INPSDIR="/pnfs/dune/persistent/users/picker24/FluxVariations"

dp_BuildUncertaintyMatrix -n ${INPSDIR}/nominal.FHC.Optimized.root -v WaterLayer,${INPSDIR}/WaterLayer_p1.FHC.Optimized.root,${INPSDIR}/WaterLayer_m1.FHC.Optimized.root -v HornCurrent,${INPSDIR}/HC_p1.FHC.Optimized.root,${INPSDIR}/HC_m1.FHC.Optimized.root -v DecayPipeR,${INPSDIR}/DecayPipeR_p1.FHC.Optimized.root,${INPSDIR}/DecayPipeR_m1.FHC.Optimized.root -E 5 -X 35 -o numu.uncerts.FHC.Optimized.root
