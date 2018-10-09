# Analysis tools usage prompts

The usage text in each executable is generated from this file, so while this
file cannot be 'out of date' relative to the usage text, it is possible for
the usage text to not reflect the current state of each executables CLI.

## `dp_RunSelection`

```
    -i <stopprocessor.root>     : TChain descriptor for input tree.
    -o <outputfile.root>        : Output file to write selected tree to.
    -v <hadr veto threshold>    : Hadronic shower veto threshold in MeV.
    -FV <fvx,y,z>               : Vertex selection fiducial volume padding
                                  inside of the non-veto active region
                                  {Default: 0,0,0}.
    -m <muon exit KE>           : Muon exit threshold KE in MeV.
    -A <EHadrVisMax_GeV>        : Acceptance cut for hadronic shower energy.
                                  Showers with more than this energy are cut
                                  and filled in by far detector MC.
    -FDproton                   : Build missing proton energy fake data. This
                                  is very hard coded, if you don't know what it
                                  is, don't use it.
```

## `dp_ProduceEfficiencyCorrector`

```
    -i <stopprocessor.root>     : TChain descriptor for input tree.
    -o <outputfile.root>        : Output file to write selected tree to.
    -v <hadr veto threshold>    : Hadronic shower veto threshold in MeV
                                  {default: 10}.
    -A <EHadrVisMax_GeV>        : Acceptance cut for hadronic shower energy.
                                  Showers with more than this energy are cut
                                  and filled in by far detector MC. Event
                                  failing this are not included in the
                                  efficiency calculation.
    -m <muon exit KE>           : Muon exit threshold KE in MeV
                                  {default: 0}.
    -b <binning descriptor>     : Energy binning descriptor that can be used
                                  to produce a perfect efficiency correction
                                  in the relevant off-axis bins.
    -FDproton                   : Build missing proton energy fake data. This
                                  is very hard coded, if you don't know what
                                  it is, don't use it.
```

## `dp_EfficiencyWeightFriendTreeBuilder`

```
    -i <stopprocessor.root>     : TChain descriptor for input tree.
                                  Should be the output of RunSelection.    
    -o <outputfile.root>        : Output file to write efficiency friend tree
                                  to.
    -E <effcorrector.root>      : File containing efficiency histograms.
    -m <muon exit KE>           : Muon exit threshold KE in MeV.
    -M <eff correction mode>    : Kinematics to use to determine the hadronic
                                  shower selection efficiency correction.
                                  {Default = 3}
                                  1: True hadronic energy, absolute off axis
                                     position.
                                  2: True hadronic energy, position within the
                                     detector.
                                  3: True non-neutron energy, position within
                                     the detector.
                                  4: Visible hadronic energy, position within
                                     the detector.
                                  5: True non-neutron energy, absolute off axis
                                     position.
```

## `dp_PRISMAnalysis`

```
    -F <fluxfitresults.root>             : Result of dp_FitFluxes to use to
                                           build observation.
    -NI <NDEvents.root>                  : TChain descriptor for input tree.
    -NF <NDFluxFile.root>                : Friend tree for -NI option containing
                                           flux throws.
    -NX <NDXSecFile.root>                : Friend tree for -NI option containing
                                           xsec throws.
    -NE <NDEffFile.root>                 : Friend tree for -NI option containing
                                           efficiency correction weights.
    -ND <NDDataFile.root>                : File containing ND data
                                           distributions.
    -NA <EHadrVis_GeV>                   : Shower acceptance cut in GeV.
    -FI <FDEvents.root>                  : TChain descriptor for FD input tree.
    -FF <FDFluxFile.root>                : Friend tree for -FI option containing
                                           flux throws.
    -FX <FDXSecFile.root>                : Friend tree for -FI option containing
                                           xsec throws.
    -FE <FDEffFile.root>                 : Friend tree for -FI option containing
                                           efficiency correction weights.
    -FD <FDDataFile.root>                : File containing FD data distribution.
    -PV <projID>                         : Kinematic projection to use
                                           {default = 3}.
                                           1: True neutrino energy
                                           2: Final state lepton energy +
                                              proton kinetic energy +
                                              non-neutron total energy
                                           3: Final state lepton energy +
                                              contained GEANT4 hadronic
                                              energy deposits (ERec).
                                           4: Final state lepton energy.
    -b <low_up:width[,low_up:width...]>  : Projection binning descriptor (GeV).
    -o <OutputFile.root>                 : Output file name.
```

## `dp_ProduceDiagnosticPlots`

```
    -i <stopprocessor.root>     : TChain descriptor for input tree.
    -o <outputfile.root>        : Output file to write selected tree to.
    -v <hadr veto threshold>    : Hadronic shower veto threshold in MeV.
    -FV <fvx,y,z>               : Vertex selection fiducial volume padding
                                  inside of the non-veto active region
                                  {Default: 0,0,0}.
    -m <muon exit KE>           : Muon exit threshold KE in MeV.
    -A <EHadrVisMax_GeV>        : Acceptance cut for hadronic shower energy.
                                  Showers with more than this energy are cut and
                                  filled in by far detector MC.
    -E <effcorrector.root>      : File containing efficiency histograms.
    -F <fluxfitresults.root>    : Result of dp_FitFluxes to use to build
                                  observation.
    -M <eff correction mode>    : Kinematics to use to determine the hadronic
                                  shower selection efficiency correction.
                                  {Default = 3}
                                  1: True hadronic energy, absolute off axis
                                     position.
                                  2: True hadronic energy, position within the
                                     detector.
                                  3: True non-neutron energy, position within
                                     the detector.
                                  4: Visible hadronic energy, position within
                                     the detector.
                                  5: True non-neutron energy, absolute off axis
                                     position.
    -b <binning descriptor>     : Energy binning descriptor that can be used to
                                  produce a perfect efficiency correction in
                                  the relevant off-axis bins.
    -S                          : Make each diagnostic plot for each detector
                                  stop and each linear combination slice.
```

## `dp_XSecVarFriendTreeBuilder`

```
    -i <stopprocessor.root>     : TChain descriptor for input tree.
                                  Should be the output of RunSelection.     
    -o <outputfile.root>        : Output file to write efficiency friend tree
                                  to.
    -C <covmatfile.root,hname>  : Input covariance matrix to throw xsec weights
                                  from.
    -N <NThrows>                : Number of throws to make.
                                  Special cases are:
                                    1: Will produce +1 sigma variation
                                    2: will produce +1 and -1 sigma variation.
    -n <MaxEvents>              : Maximum number of events to process.
    -S <Seed>                   : Used to make reproducible throws.
                                  {Default = 1}
    -D <diagonal addition>      : Covariance matrix extra diagonal error, can
                                  be used to help non-invertible matrices.
    -T <Cholesky decomp. tol.>  : Numerical tolerance on the Cholesky
                                  decomposition of the input matrix. 
```
## `dp_FluxVarFriendTreeBuilder`

```
    -i <fulldetprocess.root>             : TChain descriptor for input tree.
    -hn <weighthist.root,weighthistname> : Input file and histogram name to use
                                           for nominal flux in Enu:Lateral
                                           offset.
    -hv <weighthist.root,weighthistname> : Input file and histogram name to use
                                           for varied flux in Enu:Lateral
                                           offset.
    -o <outputfile.root>                 : Output file to write friend tree to.
```
