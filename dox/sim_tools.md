# Simulation tools usage prompts

These were correct at the time of writing: 2018-04-05

## `dp_Condenser`

```
[INFO]: Use like: dp_Condenser
	-i <inputg4arbofile>     : Input G4Py root file.
	-ir <GENIERTFile>        : Input GENIE rootracker file (index sync'd with -i argument)
	-o <outputfile>          : File to write output to
	-dmn <detxmin,ymin,zmin> : Active detector minimum (cm)
	-dmx <detxmax,ymax,zmax> : Active detector maximum (cm)
	-V <vetogap x,y,z>       : Active veto region to pad each corresponding face of the
	                           active volume with.
	                          -- N.B. If -dmn -1,-1,-1 -dmx 1,1,1 -V 0.5,0.5,0.5
	                           then the non-veto volume will be 1x1x1 cm^{3} centered on
	                           the origin.
	-P <POTPerFile>          : Adds POTPerFile information to the metadata, used for
	                           POT-normalising predicted event rates and Near/Far
	                           comparisons downstream.
	-n <NMaxEvents>          : Run no more than -n events.
	-A                       : Output all events even if they occured outside of the
	                           non-veto active volume.
	-T                       : Will add timesep branches to the output that
	                           contain all deposits ocurring more than -T <timesep>
	                           microseconds after the neutrino interaction.
	-nx <NXSteps>            : Number of x slices to break up total non-veto active
	                           region into.
	                          -- N.B. two steps will be added for the X veto gap
	                             passed to -V.
	                          -- N.B. The non-veto active region X dimension,
	                             i.e. -dmx less -dmn less times the -V should be
	                             easily divisible by this number: e.g. 10 cm steps.
	-nt <NMaxTrackSteps>     : Track final state charged lepton through up to -nt GEANT4
	                           steps.
```

## `dp_StopProcessor`

```
[USAGE]: dp_StopProcessor
	-i <input dir>                : Input directory containing condensed arbox.py output.
	-r <RunPlan.XML>              : An XML file specifying a run plan  
	                                to place stops for.
	                              -- See ${DUNEPRISMTOOLSROOT}/configs/run_plans for examples.
	-o <output.root>              : Output file name.
	-v <veto threshold MeV>       : Threshold energy deposit in veto region to pass
	                                selection {default = 10 MeV}.
	-A                            : Write out all events regardless of whether they
	                                fall within a stop.
	-n  <NEvents>                 : NMax events to write.
	-ns <NEvents>                 : N events to skip.
	-P <POTPerFile>               : POT Per input file. Overrides POT info stored in Condesed config tree.
```
