#!/bin/bash

${DUNEPRISMTOOLSROOT}/scripts/flux_scripts/FarmG4LBNE_PPFX_Makedk2nuLite.sh \
  -m ${DUNEPRISMTOOLSROOT}/configs/g4lbnf_macros/OptimizedEngineeredNov2017Review.mac \
  -p test_wppfx/v3r5p4/QGSP_BERT/OptimizedEngineeredNov2017Review/neutrino \
  --number-of-jobs 2 --number-of-protons 5000 --PPFX
