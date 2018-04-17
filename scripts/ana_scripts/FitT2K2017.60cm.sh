#!/bin/bash

dp_FitFluxes -f ~/DUNEPrismFluxesForFits/variations/nominal.FHC.uniform.root,numu_flux_2D -t osc.root,fardet_osc \
  -M 0_36:0.6 -o T2K2017.fluxfit.results.60cm.root -p -m 3 -ms 1 -of 0.01 -rg 0.1 \
  -n 10000000
