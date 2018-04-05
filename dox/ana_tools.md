# Analysis tools usage prompts

These were correct at the time of writing: 2018-04-05

## `dp_RunSelection`

```
[USAGE]: dp_RunSelection
	-i <stopprocessor.root>     : TChain descriptor for input tree.
	-o <outputfile.root>        : Output file to write selected tree to.
	-v <hadr veto threshold>    : Hadronic shower veto threshold in MeV.
	-FV <fvx,y,z>               : Vertex selection fiducial volume padding inside of the
	                              non-veto active region {Default: 0,0,0}.
	-m <muon exit KE>           : Muon exit threshold KE in MeV.
	-A <EHadrVisMax_GeV>        : Acceptance cut for hadronic shower energy. Showers with
	                              more than this energy are cut and filled in by far detector MC.
	-FDproton                   : Build missing proton energy fake data. This is very hard coded, if you don't know what it is, don't use it.
```

## `dp_ProduceEfficiencyCorrector`

```
[USAGE]: dp_ProduceEfficiencyCorrector
	-i <stopprocessor.root>     : TChain descriptor for input tree.
	-o <outputfile.root>        : Output file to write selected tree to.
	-v <hadr veto threshold>    : Hadronic shower veto threshold in MeV {default: 10}.
	-m <muon exit KE>           : Muon exit threshold KE in MeV {default: 0}.
	-b <binning descriptor>     : Energy binning descriptor that can be used to produce a perfect efficiency correction in the relevant off-axis bins.
	-FDproton                   : Build missing proton energy fake data. This is very hard coded, if you don't know what it is, don't use it.
```

## `dp_EfficiencyWeightFriendTreeBuilder`

```
[USAGE]: dp_EfficiencyWeightFriendTreeBuilder
	-i <stopprocessor.root>     : TChain descriptor for input tree.
	                              Should be the output of RunSelection.	-o <outputfile.root>        : Output file to write efficiency friend tree to.
	-E <effcorrector.root>      : File containing efficiency histograms.
	-m <muon exit KE>           : Muon exit threshold KE in MeV.
	-M <eff correction mode>    : Kinematics to use to determine the hadronic shower selection
	                              efficiency correction. {Default = 3}
	                              1: True hadronic energy, absolute off axis position.
	                              2: True hadronic energy, position within the detector.
	                              3: Visible hadronic energy, position within the detector.
```

## `dp_PRISMAnalysis`

```
[USAGE]: dp_PRISMAnalysis
	-F <fluxfitresults.root>        : Result of dp_FitFluxes to use to build observation.
	-NI <NDEvents.root>            : TChain descriptor for input tree.
	-NF <NDFluxFile.root>          : Friend tree for -NI option containing flux throws.
	-NX <NDXSecFile.root>          : Friend tree for -NI option containing xsec throws.
	-NE <NDEffFile.root>           : Friend tree for -NI option containing efficiency correction weights.
	-ND <NDDataFile.root>          : File containing ND data distributions.
	-NA <EHadrVis_GeV>             : Shower acceptance cut in GeV.
	-FI <FDEvents.root>            : TChain descriptor for FD input tree.
	-FF <FDFluxFile.root>          : Friend tree for -FI option containing flux throws.
	-FX <FDXSecFile.root>          : Friend tree for -FI option containing xsec throws.
	-FE <FDEffFile.root>           : Friend tree for -FI option containing efficiency correction weights.
	-FD <FDDataFile.root>          : File containing FD data distribution.
	-b <low_up:width[,low_up:width...]>  : ENuBinning descriptor.
	-o <OutputFile.root>                 : Output file name.
```

## `dp_ProduceDiagnosticPlots`

```
[USAGE]: dp_ProduceDiagnosticPlots
	-i <stopprocessor.root>     : TChain descriptor for input tree.
	-o <outputfile.root>        : Output file to write selected tree to.
	-v <hadr veto threshold>    : Hadronic shower veto threshold in MeV.
	-FV <fvx,y,z>               : Vertex selection fiducial volume padding inside of the
	                              non-veto active region {Default: 0,0,0}.
	-m <muon exit KE>           : Muon exit threshold KE in MeV.
	-A <EHadrVisMax_GeV>        : Acceptance cut for hadronic shower energy. Showers with
	                              more than this energy are cut and filled in by far detector MC.
	-E <effcorrector.root>      : File containing efficiency histograms.
	-F <fluxfitresults.root>        : Result of dp_FitFluxes to use to build observation.
	-M <eff correction mode>    : Kinematics to use to determine the hadronic shower selection
	                              efficiency correction. {Default = 3}
	                              1: True hadronic energy, absolute off axis position.
	                              2: True hadronic energy, position within the detector.
	                              3: Visible hadronic energy, position within the detector.
	-b <binning descriptor>     : Energy binning descriptor that can be used to produce a perfect efficiency correction in the relevant off-axis bins.
```
