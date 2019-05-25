#!/bin/bash

IDIR="finebin"

OUTPUT_FILE_NAME="DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_fine.root"

DO_PPFX="1"
DO_FOCUS="1"
DO_ALIGN="1"

for DET in "ND" "FD"; do
  for i in nu nubar; do

    if [ "${DO_PPFX}" == "1" ]; then

      dp_CombineBuiltFluxes  \
        -i "/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx/DUNEPrismFluxes/${DET}_${i}/${IDIR}/flux/Fluxes.*.root" \
        -i "/pnfs/dune/persistent/users/picker24/nominal_5E8POT_wppfx_2/DUNEPrismFluxes/${DET}_${i}/${IDIR}/flux/Fluxes.*.root" \
        --NPPFXU 100 \
        -a ${OUTPUT_FILE_NAME} \
        -D ${DET}_${i}_ppfx

    fi

    if [ "${DO_FOCUS}" == "1" ]; then
      #With focussing
      for j in p1 m1; do
        for k in WL HC DPR TargetDensity BeamSigma BeamOffsetX BeamTheta BeamThetaPhi; do

          if [ ! -e /pnfs/dune/persistent/users/picker24/Focussing/DUNEPrismFluxes/${DET}_${i}/${k}${j}/${IDIR}/flux ]; then
            continue;
          fi

        dp_CombineBuiltFluxes  \
          -i "/pnfs/dune/persistent/users/picker24/Focussing/DUNEPrismFluxes/${DET}_${i}/${k}${j}/${IDIR}/flux/Fluxes.*.root" \
          -i "/pnfs/dune/persistent/users/picker24/Focussing_2/DUNEPrismFluxes/${DET}_${i}/${k}${j}/${IDIR}/flux/Fluxes.*.root" \
          -a ${OUTPUT_FILE_NAME} \
          -D ${DET}_${i}_${k}_${j}

        done
      done
    fi

    if [ "${DO_ALIGN}" == "1" ]; then
      for j in Horn1 Horn2; do
        for k in X Y XNeg X3mm XNeg3mm; do

          if [ ! -e /pnfs/dune/persistent/users/picker24/Alignment/DUNEPrismFluxes/${DET}_${i}/${j}${k}/${IDIR}/flux ]; then
            continue;
          fi

          dp_CombineBuiltFluxes  \
            -i "/pnfs/dune/persistent/users/picker24/Alignment/DUNEPrismFluxes/${DET}_${i}/${j}${k}/${IDIR}/flux/Fluxes.*.root" \
            -i "/pnfs/dune/persistent/users/picker24/Alignment_2/DUNEPrismFluxes/${DET}_${i}/${j}${k}/${IDIR}/flux/Fluxes.*.root" \
            -a ${OUTPUT_FILE_NAME} \
            -D ${DET}_${i}_${j}_${k}Shift

        done
      done
    fi
  done
done
