#Flux tools usage prompts

The usage text in each executable is generated from this file, so while this
file cannot be 'out of date' relative to the usage text, it is possible for
the usage text to not reflect the current state of each executable's CLI.

## `dp_MakeLitedk2nu`

```
    -i /path/to/DUNE/dk2nu/files  : Can include wildcards (remember to quote
                                    to avoid shell expansion.)              
    -p /path/to/DUNE/ppfx_friend  : Adds ppfx friend tree from passed files.
                                    Can include wildcards, but it is advisable
                                    to fully specify input paths for both --i
                                    and -p options when using this, as the
                                    ordering of the friend files must be kept
                                    consistent.
    -N <Number of PPFX universes> : Defaults to 100.
    -o output.root                : File to write lite dk2nu TTree to.         
    -?                            : Display this message.

```

## `dp_BuildFluxes`

```
    -i /path/to/DUNE/dk2nu/files  : Can include wildcards (remember to quote
                                    to avoid shell expansion.)              

    -o output.root                : File to write fluxes to.                

    -m 0_1:0.5[,2,3]              : Flux window binning to calculate in mrads.

    -d 0_1:0.5[,2,3]              : Flux window binning to calculate in degrees.

    -x 0_1:0.5[,2,3]              : Flux window binning to calculate in lateral
                                    offset (m).

    -e                            : Build fluxes for specific neutrino decay
                                    parents.                                

    -b <NBins>,<Low>,<High>       : Use uniform binning for flux histograms.

    -vb 0_1:0.5[,2,3]             : Use variable binning specified by bin edges
                                    and step ranges.                        

    -h <Height=0>                 : Height of flux plane (cm).

    -n <NMaxNeutrinos>            : Only loop over -n neutrinos.    

    -z <ZDist>                    : Z distance of flux plane from target (cm).

    -P                            : Only use each decaying parent once.       
    -S <species PDG>              : Only fill information for neutrinos of a
                                    given species.
    -L                            : Expect dk2nulite inputs.
    --PPFX                        : Expects to be able to read PPFX weight
                                    branches from input dk2nu files (note, only
                                    dk2nulite files will likely have these
                                    branches.)
    --NPPFXU <int>                : The number of PPFX universes, defaults to
                                    100.
    -?                            : Display this message.
```

## `dp_CombineBuiltFluxes`

```
    -i <Input search pattern>  : Search pattern to find input files. Can be
                                 specified multiple times.
    -o <Output file name>      : File to write combined output to.
    --NPPFXU <int>             : The number of PPFX universes. Defaults to 0.
    -?                         : Display this message.
```

## `dp_FitFluxes`

```
  Input options:                                                          

    -f <ROOT file,FluxHist2DName>      : Input 2D flux histogram, Y bins
                                         correspond to different fluxes.

    -MX  <nbins to merge>              : Merge neutrino energy bins before
                                         splitting into fluxes.

    -M  <OA1>:<OA_W>,<OA2>_<OAN>:<OA_W>,<OAN+1>:<OA_W>,...
                                       : Merge bins in off axis flux positions.
                                         Each position or position range
                                         specifies a slice width. The
                                         corresponding absolute slice ranges
                                         must match up to merge-able Y bin
                                         edges from histogram passed to -f.
                                         You will be notified if they don't.

    -A <FitOutput.root[,dirname]>      : Start a new fit from the results of an
                                         old fit. Optional dirname corresponds
                                         to -d option. (Tip: set -n 0 to apply
                                         previous results to new inputs without
                                         running a fit.)

  Output options:                                                         
    -[o|a] <ROOT file>                 : The output root file. Using -o will
                                         overwrite a file of the same name,
                                         -a will append the fit result to   
                                         the file.                          

    -d <directory name>                : If passed, fit result will be put  
                                         into a subdirectory of the root    
                                         file.                              

  Target options:                                                         
    -t <ROOT file,hist name>           : The histogram of the target flux to
                                         fit to. This file should contain an
                                         oscillation parameter config tree
                                         generated by dp_OscillateFlux.

    -g <mean,width>                    : Use a gaussian target distribution
                                         instead of a target flux shape.    

  Fitter options:                                                       
    -n <MaxCalls=50000>                : The maximum number of Likelihood       
                                         evaluations before giving up the   
                                         fit.

    -c <CoeffLimit=30>                 : Parameter limits of flux component
                                         coefficients.                      

    -T <Tolerance=1E-5>                : Minuit2 EDM tolerance.
  Figure of merit options:                                                
    -C                                 : Use NuPrism tools Chi2.

    -rg <regularisation factor>        : Adds neighbouring coefficient
                                         regularisation.                    

    -l <min val>,<max val>             : Fit between min and max. Outside
                                         of this range, -m determines
                                         behaviour.             

    -p                                 : Fit between the first and third
                                         oscillation peaks (Only useful with -t)

    -m <0,1,2,3>                       : Out of range behaviour.            
                                         0: Ignore out of range bins.      
                                         1: Force out of range bins to 0.  
                                         2: Exponential decay outside fit  
                                            region. Decay rate is          
                                            determined by -ed.
                                         3: Gaussian decay outside fit  
                                            region. Decay width is          
                                            determined by -ed.

    -ed <decay rate>                   : For -m [2,3], controls decay rate.    
                                         Default = 3, larger is faster     
                                         decay.

    -of <out of range factor>          : Allow out of range to contribute less
                                         to the FOM by this factor.

    -ms <out of range side>            : Controls which 'side' of the fit
                                         regions are controlled by -m
                                         0: Include both low and high out of
                                            range E,
                                         1: Include low E
                                         2: Include high E.
```

## `dp_OscillateFlux`

```
    -p <p1,p2,p3,p4,p5,p6>              : A comma separated list of oscillation
                                          parameters:
                                           sin2(theta12) {default = 0.825  }
                                           sin2(theta13) {default = 0.10   }
                                           sin2(theta23) {default = 1.0    }
                                           dm12          {default = 7.9e-5 }
                                           dm23          {default = 2.5e-3 }
                                           dcp           {default = 0.0    }
    -d <dipangle>                       : Beam dip angle in degrees. Used to
                                          calculate the oscillation baseline.
                                          N.B. currently the far point is
                                          assumed to be on the surface.
    -i <file.root,histname>             : Input neutrino flux histogram to apply
                                          oscillation weights. X-axis assumed
                                          to be E_{nu} in GeV.
    -n <osc from, osc to>               : PDG MC particle codes describing
                                          neutrino oscillation channel (e.g.
                                          14,14 for muon neutrino
                                          disappearance).
    -L <baseline>                       : Alternative to -d, specify baseline in
                                          km. Used to determine the beam dip
                                          angle needed to give this baseline
                                          (and thus get the matter effect
                                          approximately correct through Prob3++
                                          use of earth-like matter
                                          distribution.)
  Output options:                                                         
    -[o|a] <ROOT file>                  : The output root file. Using -o will
                                          overwrite a file of the same name,
                                          -a will append the fit result to   
                                          the file.                          

    -D <directory name>                 : If passed, fit result will be put  
                                          into a subdirectory of the root    
                                          file.      
```

## `dp_FluxSmusher`

```
  Input options:   
    -i <ROOT file,FluxHist1DName>  : Input histogram, multiple instances can be
                                     supplied
    -f <ROOT file,FluxHist2DName>  : Input 2D flux histogram, Y bins
                                     correspond to different fluxes.
    -M  <OA1>:<OA_W>,<OA2>_<OAN>:<OA_W>,<OAN+1>:<OA_W>,...
                                   : Merge bins in off axis flux positions.
                                     Each position or position range
                                     specifies a slice width. The
                                     corresponding absolute slice ranges
                                     must match up to merge-able Y bin
                                     edges from histogram passed to -f.
                                     You will be notified if they don't.
  Output options:                                                         
    -[o|a] <ROOT file>             : The output root file. Using -o will
                                     overwrite a file of the same name,
                                     -a will append the fit result to   
                                     the file.                          
```

## `dp_BuildUncertaintyMatrix`

```
    --fhicl                       : FHiCL configuration file, see
                                    ${DUNEPRISMTOOLSROOT}/fcl/flux_uncertainty_unputs.fcl
                                    for an example.         
    -?                            : Display this message.

```
