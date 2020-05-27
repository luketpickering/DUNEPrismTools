#!/bin/bash

IDIR="uncertbin_offaxis"

OUTPUT_FILE_NAME="DUNE_Flux_OffAxis_Nov2017Review_syst_shifts_uncertbin_offaxis_280kAOnAxis.root"

DO_PPFX="1"
DO_PPFX_COMPONENT_VARIATIONS="0"
DO_FOCUS="1"
DO_ALIGN="1"
DO_HC280KA="1"
HC280KA_DIR="280KA"
if [ "${DO_PPFX}" == "1" ]; then
  HC280KA_DIR="280KA_wppfx"
fi

for DET in "ND" "FD"; do
  for i in nu nubar; do

    if [ "${DO_PPFX}" == "1" ]; then

      dp_CombineBuiltFluxes  \
        -i "/pnfs/dune/persistent/users/${USER}/nominal_5E8POT_wppfx/DUNEPrismFluxes/${DET}_${i}/${IDIR}/flux/Fluxes.*.root" \
        --NPPFXU 100 \
        -a ${OUTPUT_FILE_NAME} \
        -D ${DET}_${i}_ppfx

    fi

    if [ "${DO_PPFX_COMPONENT_VARIATIONS}" == "1" ]; then

      dp_CombineBuiltFluxes  \
        -i "/pnfs/dune/persistent/users/${USER}/nominal_2.5E8POT_wallppfx/DUNEPrismFluxes/${DET}_${i}/${IDIR}/flux/Fluxes.*.root" \
        --NPPFXU 100 --ReadPPFXAllWeights \
        -a ${OUTPUT_FILE_NAME} \
        -D ${DET}_${i}_ppfx_allw

    fi

    if [ "${DO_FOCUS}" == "1" ]; then
      #With focussing
      for j in p1 m1; do
        for k in WL HC DPR TargetDensity BeamSigma BeamOffsetX BeamTheta BeamThetaPhi; do

          if [ ! -e /pnfs/dune/persistent/users/${USER}/Focussing/DUNEPrismFluxes/${DET}_${i}/${k}${j}/${IDIR}/flux ]; then
            continue;
          fi

        dp_CombineBuiltFluxes  \
          -i "/pnfs/dune/persistent/users/${USER}/Focussing/DUNEPrismFluxes/${DET}_${i}/${k}${j}/${IDIR}/flux/Fluxes.*.root" \
          -a ${OUTPUT_FILE_NAME} \
          -D ${DET}_${i}_${k}_${j}

        done
      done
    fi

    if [ "${DO_ALIGN}" == "1" ]; then
      for j in Horn1 Horn2; do
        for k in X Y XNeg X3mm XNeg3mm; do

          if [ ! -e /pnfs/dune/persistent/users/${USER}/Alignment/DUNEPrismFluxes/${DET}_${i}/${j}${k}/${IDIR}/flux ]; then
            continue;
          fi

          dp_CombineBuiltFluxes  \
            -i "/pnfs/dune/persistent/users/${USER}/Alignment/DUNEPrismFluxes/${DET}_${i}/${j}${k}/${IDIR}/flux/Fluxes.*.root" \
            -a ${OUTPUT_FILE_NAME} \
            -D ${DET}_${i}_${j}_${k}Shift

        done
      done
    fi

    if [ "${HC280KA_DIR}" == "1" ]; then
      if [ ! -e /pnfs/dune/persistent/users/${USER}/${HC280KA_DIR}/DUNEPrismFluxes/${DET}_${i}/${IDIR}/flux ]; then
        echo "Cannot find /pnfs/dune/persistent/users/${USER}/${HC280KA_DIR}/DUNEPrismFluxes/${DET}_${i}/${IDIR}/flux"
        continue;
      fi

      dp_CombineBuiltFluxes  \
        -i "/pnfs/dune/persistent/users/${USER}/${HC280KA_DIR}/DUNEPrismFluxes/${DET}_${i}/${IDIR}/flux/Fluxes.*.root" \
        -a ${OUTPUT_FILE_NAME} \
        -D ${DET}_${i}_HC_${CURR} --NPPFXU 100

    fi
  done
done
