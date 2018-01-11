# NuPRISM / E61 scripts hacked for DunePRISM
Very messy, lots of hordcoded paths, etc

 - fit_gaussian.cc : Looks for a 'duneprism_spectra.root' produced by makeDune2Dhist.C finds the linear combination that fits a Gaussian with defined mean and width.
 - fit_spectrum.cc : Fits fluxes in 'duneprism_spectra.root' to oscillated flux histogram 'numu_fluxosc_forplots' in file pointed at by 'fDuneFlux'
 - makeDune2Dhist.C : Produces 'duneprism_spectra.root' file from off-axis fluxes produced by the tools in 'offAxis'

 - offAxis : Tools to produce oscillated and unoscillated fluxes at different detector positions. Needs beam simulation files.
