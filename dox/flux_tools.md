# Flux tools usage prompts

These were correct at the time of writing: 2018-04-05

## `dp_BuildFluxes`

```
[USAGE]: dp_BuildFluxes
	-i /path/to/DUNE/dk2nu/files  : Can include wildcards (remeber to quote
	                                to avoid shell expansion.)              

	-o output.root                : File to write fluxes to.                

	-m 0_1:0.5[,2,3]              : Flux window binning to calculate in mrads.

	-d 0_1:0.5[,2,3]              : Flux window binning to calculate in degrees.

	-x 0_1:0.5[,2,3]              : Flux window binning to calculate in lateral offset (m).

	-e                            : Build fluxes for specific neutrino decay
	                                parents.                                

	-b <NBins>,<Low>,<High>       : Use uniform binning for flux histograms.

	-vb 0_1:0.5[,2,3]             : Use variable binning specified by bin edges
	                                and step ranges.                        

	-h <Height=0>                 : Height of flux plane (cm).

	-n <NMaxNeutrinos>            : Only loop over -n nus.    

	-z <ZDist>                    : Z distance of flux plane from target (cm).

	-P                            : Only use each decaying parent once.       
	-S <species PDG>              : Only fill information for neutrinos of a given species.
	-L                            : Expect dk2nulite inputs.
```

## `dp_CombineBuiltFluxes`

```
[USAGE]: dp_CombineBuiltFluxes
	-i <Input search pattern>  : Search pattern to find input files. Can be specified multiple times.

	-o <Output file name>      : File to write combined output to.
	-?                         : Display this message.
```

## `dp_FitFluxes`

```
[USAGE]: dp_FitFluxes
  Input options:                                                          

	-f <ROOT file,FluxHist2DName>      : Input 2D flux histogram, Y bins
	                                     correspond to different fluxes.

	-MX  <nbins to merge>              : Merge neutrino energy bins before splitting into
	                                     fluxes.

	-M  <OA1>:<OA_W>,<OA2>_<OAN>:<OA_W>,<OAN+1>:<OA_W>,...
 	                                   : Merge bins in off axis flux positions. Each position
	                                     or position range specifies a slice width. The
	                                     corresponding absolute slice ranges must match up to
	                                     merge-able Y bin edges from histogram passed to -f.
	                                     You will be notified if they don't.

	-A <FitOutput.root[,dirname]>      : Start a new fit from the results of an
	                                     old fit. Optional dirname corresponds to
	                                     -d option. (Tip: set -n 0 to apply previous
	                                     results to new inputs without running a fit.)

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
	                                     oscillation parameter config tree generated
	                                     by dp_OscillateFlux.

	-g <mean,width>                    : Use a gaussian target distribution
	                                     instead of a target flux shape.    

  Fitter options:                                                       
	-n <MaxCalls=50000>                : The maximum number of MINUIT       
	                                     evaluations before giving up the   
	                                     fit.

	-c <CoeffLimit=30>                 : Parameter limits of flux component
	                                     coefficients.                      

  Figure of merit options:                                                
	-C                                 : Use NuPrism tools Chi2.

	-rg <regularisation factor>        : Adds neighbouring       coefficient
	                                     regularisation.                    

	-l <min val>,<max val>             : Fit between min and max. Outside
	                                     of this range, -mdetermines
	                                     behavior.             

	-p                                 : Fit between the firstand third oscillation
	                                     peaks (Only useful with -t)

	-m <0,1,2,3>                       : Out of range behavior.            
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

	-of <out of range factor>          : Allow out of range to contribute less to the FOM
 	                                     by this factor.

	-ms <out of range side>            : 0 = Include both low and high out of range E,
	                                     1 = include low E, 2 = include high E.
```
