#!/bin/bash

# dp_ProduceDiagnosticPlots \
#   -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30/*.root" \
#   -v 20 \
#   -o DiagnosticPlots.4mwide.0cmFV.EHadrTrueAbsPos.NoAcceptance.root \
#   -E EffCorrector.4mwide.0cmFV.NoAcceptance.root \
#   -F FluxFit.root \
#   -b 0_10:0.1 \
#   -M 1 -S
#
# dp_ProduceDiagnosticPlots \
#   -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC/Processed.2018-03-30/*.root" \
#   -v 20 \
#   -o DiagnosticPlots.4mwide.0cmFV.ENoneNeutronHadr.NoAcceptance.root \
#   -E EffCorrector.4mwide.0cmFV.NoAcceptance.root \
#   -F FluxFit.root \
#   -b 0_10:0.1 \
#   -M 3 -S
#
# dp_ProduceDiagnosticPlots \
#   -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07/*.root" \
#   -v 20 \
#   -o DiagnosticPlots.7mwide.100cmFV.ENoneNeutronHadr.NoAcceptance.root \
#   -E EffCorrector.7mwide.100cmFV.NoAcceptance.root \
#   -F FluxFit.root \
#   -b 0_10:0.1 \
#   -M 3 -FV 100,0,0 -S

dp_ProduceDiagnosticPlots \
  -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07/*.root" \
  -v 20 \
  -o DiagnosticPlots.7mwide.100cmFV.ENoneNeutronHadr.detpos_corr.NoAcceptance.root \
  -E EffCorrector.7mwide.100cmFV.NoAcceptance.root \
  -F FluxFit.root \
  -b 0_10:0.1 \
  -M 3 -FV 100,0,0 -S

dp_ProduceDiagnosticPlots \
  -i "/mnt/research/NuInt/Dune_Flux/OptimizedEngineeredSept2017Review/Analysis/FHC_overlaps7m/Processed.2018-04-07/*.root" \
  -v 20 \
  -o DiagnosticPlots.7mwide.100cmFV.ENoneNeutronHadr.abspos_corr.NoAcceptance.root \
  -E EffCorrector.7mwide.100cmFV.NoAcceptance.root \
  -F FluxFit.root \
  -b 0_10:0.1 \
  -M 5 -FV 100,0,0 -S
