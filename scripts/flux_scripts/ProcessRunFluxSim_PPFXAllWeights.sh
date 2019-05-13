#!/bin/bash

# for i in nu nubar; do
for i in nu; do

    if [ -e /pnfs/dune/persistent/users/picker24/nominal_2.5E8POT_wallppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i} ]; then
      echo "[INFO]: Not regenerating /pnfs/dune/persistent/users/picker24/nominal_2.5E8POT_wallppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i}"
      continue
    fi

    if [ ! -e ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode.mac ]; then
      echo "[INFO]: No such g4lbnf macro OptimizedEngineeredNov2017Review_${i}mode.mac Not running."
      continue
    fi

    ${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
      -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review_${i}mode.mac \
      -p nominal_2.5E8POT_wallppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/${i} \
      --number-of-jobs 1 --PPFX --Use-All-PPFX-Weights \
      --expected-disk 2GB \
      --expected-mem 1.9GB \
      --expected-walltime 4h

done
