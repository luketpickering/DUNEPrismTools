# DUNEPrismTools

## Overview

This package contains tools for producing inputs for, and developing a
PRISM-style analysis. These are split into three distinct stages: neutrino
flux input generation, neutrino interaction simulation and final state
propagation, and event selection, systematic propagation, and PRISM prediction
production.

The structure of the package favors using multiple, focussed tools that perform
specific tasks and pass meta-data about inputs, processing, and outputs to
subsequent processing via small `ROOT` `TTree`s. Each tool should provide
a detailed description of its usage if interrogated with `-?`. See the `app`
directory documentation linked in each section below to see an example of
the usage text for each app.

Precise terminology used in the analysis is highlighted in the below
description.

### PRISM analysis

The PRISM analysis technique uses linear combinations of observables taken
under exposure to different neutrino fluxes to build a prediction of an
observable taken under some chosen neutrino flux. The most directly useful
example uses linear combinations of measurements taken at different off-axis
angles in the near detector of an LBL experiment to predict distributions of
far detector observables under different neutrino oscillation hypotheses.
For more information on PRISM-style analyses see: [arxiv:1412.3086](https://arxiv.org/pdf/1412.3086.pdf).

### Neutrino flux tools

#### Overview

This portion of the tool chain provides tools for interpreting and reducing the
output of the beam line simulation,
[`g4lbnf`](https://cdcvs.fnal.gov/redmine/projects/lbne-beamsim/wiki), as well
as performing the flux fits, which are central to the PRISM analysis technique.

App usage text [here](flux_tools.md).

#### Beam line simulation

The flux prediction for [LBNF](https://lbnf.fnal.gov/) is provided by the
`g4lbnf` package. To help interpreting the results of the beam line simulation,
two tools are provided.

  * `dp_BuildFluxes`: Reads [`dk2nu`](https://cdcvs.fnal.gov/redmine/projects/dk2nu/wiki)
  output from `g4lbnf` and produces two dimensional flux predictions for a given
  set of flux window positions.
  * `dp_CombineBuiltFluxes`: Combines the results of multiple `dp_BuildFluxes`
  runs. Useful for sharing the load of processing a large number of `dk2nu`
  files.

To propagate the results of the flux uncertainties, a flux covariance matrix
helper was developed. However, it is not build by default, it is now out of
date and needs re-working, the implementation is defined in
`app/flux_tools/BuildUncertaintyMatrix.cxx`.

#### `dk2nu` tools

The raw `d2knu` files that are produced by `g4lbnf` contain a large amount of
information useful for beam line systematic uncertainty propagation. While this
information is needed when assessing the effect of hadron production errors,
it is not needed for a large part of the input generation used for the PRISM
analysis. The `dp_MakeLitedk2nu` tool strips out all the information from
`dk2nu` files that is not needed for `dp_BuildFluxes`. This results in a factor
of ten file size reduction.

#### PRISM flux fitting

Two additional tools are provided for the use of these flux inputs in the
PRISM analysis.

  * `dp_OscillateFlux`: Applies oscillation weights to an input neutrino flux
  given a set of PMNS oscillation parameters, a beam dip angle, and an
  oscillation channel.
  * `dp_FitFluxes`: Fits two dimensional flux predictions to some input target
  neutrino flux (usually produced by `dp_OscillateFlux` or some target
  gaussian flux).

The result of this toolchain is a set of flux window definitions and associated
linear combination coefficients that can be used to define and weight near
detector measurements to build a observable prediction under the target
flux input to `dp_FitFluxes`.

**Slice** : The term 'slice' is used consistently throughout to refer to a flux
window that defines a 'measurement' used in a linear combination of
measurements. *N.B.* This is as opposed to a **stop**, which will be introduced
later.

### Neutrino interaction and final state propagation simulation tools

#### Overview

The main tool used in this part of the analysis chain is external to this
package. As a result, some of the most useful part of this package related to
simulation are configuration `XML` and `GDML` files, and shell scripts that
can be used to submit and execute GENIE simulation jobs on the grid
(specifically tailored to submission from FNAL front end machines).

To reduce detector simulation development time and re-processing, the simulation
currently generates neutrino interactions on a very wide block of liquid argon
(LAr) and a post-processing stage splits the output of this simulation up into
more realistically sized detector observations taken at a range of off-axis
angles.

App usage text [here](sim_tools.md).

#### Interaction simulation

The initial GENIE run takes `d2knu` input beam line simulation files as
introduced in the 'Neutrino flux tools' section. It takes a flux window
definition (*e.g.* `configs/flux_windows/DUNEPrismFluxWindow.xml`) that also
defines a translation of beam line simulation coordinates to flux window
coordinates. It takes a detector geometry definition in the `GDML` file format
(*e.g.* Near detector: `configs/gdml/DUNEPrismLArBox.geom.manual.gdml`,
  Far detector: `configs/gdml/OnAxisLArBox_FD_25.8mx22.6mx56.6m.gdml`).

#### Final state propagation

**N.B.** This step is not currently runnable using tools in this package.
This should be rectified, this step should be replaced by a similar step that
uses the [`edep-sim`](https://github.com/ClarkMcGrew/edep-sim) tool set.

The output of GENIE is converted from the GHEP format to the `numi_rootracker`
event format and passed to a python script, [`argon box`](https://github.com/calcuttj/argon_box)
which propagates all final state particles through a large cuboid of LAr using a
GEANT4 simulation.

#### Energy deposit summary

The output of the GEANT4 propagation is large for even a few million events.
It was decided that a reduced, intermediate format should be used that allows
detector stops to be arbitrarily positioned after the fact, but didn't keep all
of the particle tracking information output by GEANT4.

**Stop** : The term 'stop' is used consistently throughout to describe the
active extent and veto region of a detector at a given off-axis position.
It is characterized by an absolute minimum corner position and an absolute
maximum corner position for the active extent and a veto region extent in each
dimension to be removed from each of the six outer faces of the active cuboid
of LAr. *N.B.* This is as opposed to a 'slice' as defined above.

**Veto region** : The veto region of a detector stop is used to select neutrino
interactions that have a well-sampled hadronic shower. Interactions that leave
a large amount of hadronic energy near the edge of the detector were likely not
well contained and thus the total available hadronic energy not well sampled.
The total energy deposited within this volume by final state particles that are
not the primary final state lepton or descendent particles thereof is used
downstream in the analysis to select 'well reconstructed events'. *N.B.*
While interactions occurring within the veto region are not selected, the
vertex selection 'fiducial volume' can be smaller than the non-veto active
extent of a detector stop.

`dp_Condenser` reads the full GEANT4 input and outputs a binned summary of
energy deposits left in the LAr. The exact binning is configurable at runtime,
but by default a 40 m wide block of LAr is used in the interaction simulation,
and each detector stop by default uses a 50 x 50 x 50 cm veto region. The Y and
Z (vertical and ~beam direction) dimensions are separated into 3 bins: low
dimension veto region, dimension non-veto region, high dimension veto region.
The X (off-axis direction) dimension is separated into many bins: most off axis
stop low X veto region, N bins allowing arbitrary stop and fiducial
volume placement, on axis stop high X veto region. By default, the 40 m wide
block of LAr has 50 cm removed from each side for the extremal veto regions and
the remaining 39 m of active volume is separated in 390 bins to allow detector
stop and vertex selection placement as granularly as 10 cm.

In each of these bins, the deposits energy deposits simulated as left by
final-state, nuclear-leaving particles are summed over. The deposits are
separated by particle class (proton, neutron, charged pion, neutral pion, other)
and by whether the energy deposit was left by the original particle from GENIE
or by some descendent particle. To facilitate studies of timing, it is also
possible to choose a time, after which, energy deposits are separated into
a different set of positional bins. An example might be to use a predicted
drift window in the LAr detector to try and capture how much energy is missed
due to late deposits (such as neutron thermalising).

To facilitate studies of exiting muon tagging, primary final state muons are
fully tracked through the liquid argon so that their exiting position and
momentum from a chosen detector stop can be determined post-condensing.

By default, interactions occurring outside of any possible non-veto region are
ignored by `dp_Condenser`.

The information contained within the output of this stage is declared in
`src/persistency/CondensedDepositsTreeReader.hxx`.

#### Detector stop positioning

`dp_StopProcessor` reads the condensed output, places detector stops and
integrates over the veto and non-veto region deposit bins to provide an event
summary that can be used to perform the downstream PRISM analysis.
The information output by the stop processor is well documented in
`src/persistency/DepositsSummaryTreeReader.hxx`.

The stop definitions, or 'run plans', are defined by a simple `xml` format.
An example of a run plan where the non-veto active regions are used to
contiguously sample from 2 m to -37.5 m off axis can be seen in
`configs/run_plans/RunPlan.39mLAr.4mx3mx5mActive.xml`.

It is likely that a fiducial volume definition will be used to restrict the
vertex selection region to one of high hadronic shower selection efficiency,
this reduces the magnitude of corrections made in the PRISM analysis.
*N.B.* The fiducial volume definition is not handled at this step, but the
choice of overlapping or non-overlapping non-veto regions allows different
choices of fiducial volumes downstream.

In the case of a restricted vertex selection fiducial volume, to get contiguous
selection regions over the full range of  absolute  off-axis positions,
the non-veto regions will necessarily overlap. An example of a run plan where
the non-veto active regions overlap can be seen in
`configs/run_plans/RunPlan.39mLAr.4mx3mx5mActive_overlaps.xml`.

When the non-veto active regions of two or more detector stops overlap, an
interaction occurring in an overlap volume must be placed in one of the stops
for subsequent analysis. The choice of stop is weighted by the POTExposure
attribute defined in the run plan configuration files. The event is then
given a 'stop weight' which accounts for the fact that in a full simulation,
or a real data run, both stops would have 'seen' an equivalent event. Usually
these overlaps are by two stops with equivalent POT exposure, so an event will
get randomly placed in either stop and recieve a 'stop weight' of 2.

The output of this stage contains the `DepositsSummaryTree` which is the input
format for the entire analysis chain.

### PRISM analysis tools

#### Overview

The analysis tools are split into four stages: selection, efficiency correction,
Systematic assessment, and 'analysis'.

App usage text [here](ana_tools.md).

#### Selection

The `dp_RunSelection` tool
