#ifndef DEPOSITSUMMARYTREEREADER_HXX_SEEN
#define DEPOSITSUMMARYTREEREADER_HXX_SEEN

#include "TChain.h"
#include "TObjString.h"

#include <string>

/// Energy deposit and GENIE passthrough output tree
struct DepositsSummary {
  DepositsSummary();
  DepositsSummary(std::string const &treeName, std::string const &inputFiles,
       double timesep_us = 0xdeadbeef, bool IsLite = false);

  TChain *tree;
  UInt_t NFiles;
  UInt_t NEntries;
  UInt_t CEnt;

  void GetEntry(UInt_t e);

  UInt_t GetEntry();
  UInt_t GetEntries();

  double timesep_us;
  bool IsLite;

  ~DepositsSummary();

  size_t GetNPassthroughParts();

  std::pair<Int_t, Double_t *> GetPassthroughPart(size_t i);

  bool AddPassthroughPart(Int_t PDG, Double_t *fourmom);

  static const Int_t kNMaxPassthroughParts = 100;

  ///\brief The detector stop number used, refer to input xml for stop offsets.
  ///
  /// N.B. When overlapping stops are defined the event is randomly placed
  /// within one of the overlapping stops at the interaction position. The
  /// choice is weighted by the POTExposure branch in the input run plan
  /// xml.
  int stop;

  ///\brief Event weight used to account for overlapping stops only giving each
  /// event to a single stop.
  ///
  /// If not using a runplan with overlapping stops, this should always be 1.
  double stop_weight;
  ///\brief Event weight used to account for the POT exposure specified in the
  /// runplan.
  ///
  /// Applying this weight will rescale event from generated POT to correct POT
  /// given the runplan.
  double POT_weight;

  /// [GENIE P/T]:  The vertex 3-position in cm
  double vtx[3];
  ///\brief [GENIE P/T]:  The X position of the vertex relative to the centre of
  /// a
  /// stop in cm.
  double vtxInDetX;
  /// [GENIE P/T]:  The X offset of stop, stop in cm.
  double XOffset;
  /// [GENIE P/T]:  The GENIE interaction code (Full interactions tring).
  TObjString *EventCode;
  ///\brief [GENIE P/T]:  The GENIE interaction code (integer).
  ///
  /// * 1 : QE
  /// * 2 : MEC/2p2h
  /// * 3 : RES
  /// * 4 : DIS
  /// * 5 : COH
  /// * 6 : nu-e elastic
  /// * 7 : IMD
  Int_t GENIEInteractionTopology;
  /// [GENIE P/T]:  The 4-momentum of the incident neutrino in detector
  /// coordinates.
  double nu_4mom[4];
  /// [GENIE P/T]:  The elasticity of the interaction.
  double y_True;
  /// [GENIE P/T]:  The square 4-momentum transfer of the interaction.
  double Q2_True;
  /// [GENIE P/T]:  The full 4-momentum transfer of the interaction.
  double FourMomTransfer_True[4];
  ///\brief [GENIE P/T]:  The reconstructed invariant mass
  ///
  /// N.B. This assumes that the target nucleon was at rest and will not be the
  /// same W as thrown by the event generator during the cross-section
  /// calculation (single-pion production).
  double W_Rest;

  ///\brief [INTERNAL]:  The number of final state particles from the
  /// generator.
  ///
  /// N.B. It is very unlikely that a user should use this member, it is for
  /// the TTree to be able to store variable length arrays internally. To access
  /// the particle stack, use the GetNPassthroughParts and GetPassthroughPart
  /// member functions.
  Int_t NFSParts;
  ////brief [INTERNAL]:  The PDG MC codes of all final state particles.
  ///
  /// N.B. It is very unlikely that a user should use this member, it is for
  /// the TTree to be able to store variable length arrays internally. To access
  /// the particle stack, use the GetNPassthroughParts and GetPassthroughPart
  /// member functions.
  Int_t FSPart_PDG[kNMaxPassthroughParts];
  ///\brief [INTERNAL]:   The number of entries in the four momentum array.
  ///
  /// N.B. It is very unlikely that a user should use this member, it is for
  /// the TTree to be able to store variable length arrays internally. To access
  /// the particle stack, use the GetNPassthroughParts and GetPassthroughPart
  /// member functions.
  Int_t NFSPart4MomEntries;
  ///\brief [INTERNAL]:  A flattened array of the final state particle
  /// 4-momenta.
  ///
  /// N.B. It is very unlikely that a user should use this member, it is for
  /// the TTree to be able to store variable length arrays internally. To access
  /// the particle stack, use the GetNPassthroughParts and GetPassthroughPart
  /// member functions.
  Double_t FSPart_4Mom[kNMaxPassthroughParts * 4];

  /// [GENIE P/T]:  The PDG MC code of the neutrino.
  int nu_PDG;
  ///\brief [GENIE P/T]:  The PDG MC code of the primary lepton
  ///
  /// i.e. the one that was born when the neutrino shed/absorbed an exchange
  /// boson
  int PrimaryLepPDG;
  /// [GENIE P/T]:  The 4-momentum of the primary lepton
  double PrimaryLep_4mom[4];

  /// [GENIE P/T]:  The number of final state leptons in the event.
  int NLep;
  /// [GENIE P/T]:  The number of final state neutral pions in the event.
  int NPi0;
  /// [GENIE P/T]:  The number of final state charged pions in the event.
  int NPiC;
  /// [GENIE P/T]:  The number of final state protons in the event.
  int NProton;
  /// [GENIE P/T]:  The number of final state neutrons in the event.
  int NNeutron;
  /// [GENIE P/T]:  The number of final state photons in the event.
  int NGamma;
  ///\brief [GENIE P/T]:  The number of final state particles in the
  /// event with 1000 < abs(PDG) < 9999
  ///
  /// N.B. These correspond to baryonic resonance particles and arose from a
  /// proton or neutron. N times the nucleon mass should probably be removed
  /// from these events.
  int NBaryonicRes;
  ///\brief [GENIE P/T]:  The number of anti-matter nucleons in the event.
  ///
  /// Note that the anti-nucleons will will also be counted by NNeutron and
  /// NProton.
  int NAntiNucleons;
  ///\brief [GENIE P/T]:  The number of final state other particles in the
  /// event.
  ///
  /// N.B. These do not include GENIE bindinos or nuclear PDG codes.
  /// By eye, these are most often Kaons.
  int NOther;

  ///\brief [GENIE P/T]: The total kinetic energy of all neutral pions at the
  /// end of the GENIE simulation.
  double EKinPi0_True;
  ///\brief [GENIE P/T]: The total mass energy of all neutral pions at the
  /// end of the GENIE simulation.
  double EMassPi0_True;
  ///\brief [GENIE P/T]: The total kinetic energy of all charged pions at the
  /// end of the GENIE simulation.
  double EKinPiC_True;
  ///\brief [GENIE P/T]: The total mass energy of all charged pions at the
  /// end of the GENIE simulation.
  double EMassPiC_True;
  ///\brief [GENIE P/T]: The total kinetic energy of all protons at the
  /// end of the GENIE simulation.
  double EKinProton_True;
  ///\brief [GENIE P/T]: The total mass energy of all protons at the
  /// end of the GENIE simulation.
  ///
  /// N.B. It is most often the case that the mass energy of nucleons was not
  /// created during the neutrino interaction or subsequent cascade. A proxy
  /// reconstructed neutrino energy will often not use this energy.
  double EMassProton_True;
  ///\brief [GENIE P/T]: The total kinetic energy of all neutrons at the
  /// end of the GENIE simulation.
  double EKinNeutron_True;
  ///\brief [GENIE P/T]: The total mass energy of all neutrons at the
  /// end of the GENIE simulation.
  ///
  /// N.B. It is most often the case that the mass energy of nucleons was not
  /// created during the neutrino interaction or subsequent cascade. A proxy
  /// reconstructed neutrino energy will often not use this energy.
  double EMassNeutron_True;
  ///\brief [GENIE P/T]: The total energy of all photons at the
  /// end of the GENIE simulation.
  double EGamma_True;
  ///\brief [GENIE P/T]: The total energy of all other particles at the
  /// end of the GENIE simulation.
  ///
  /// N.B. These do not include GENIE bindinos or nuclear PDG codes.
  /// By eye, these are most often Kaons or Lambdas.
  double EOther_True;
  ///\brief [GENIE P/T]: The total energy of all non-primary-leptons at the
  /// end of the GENIE simulation. (i.e. hadrons)
  double Total_ENonPrimaryLep_True;
  ///\brief [GENIE P/T]: The KE of all nucleons and total energy of all other
  /// non-primary-leptons at the end of the GENIE simulation. (i.e. hadrons)
  double ENonPrimaryLep_KinNucleonTotalOther_True;

  ///\brief [GENIE P/T]: The summed three momentum of all final state particles
  /// from GENIE.
  double TotalFS_3mom[3];

  ///\brief [GENIE P/T]: The KE of all protons and total energy of all other
  /// particles (including primary proton), with 938 MeV removed for all
  /// baryonic resonances found.
  double ERecProxy_True;

  ///\brief [GEANT4]: The total 'early' lepton energy deposited within the stops
  /// fiducial volume.
  ///
  /// N.B. Unlike other branches, this does *not* perform descendent roll up.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double LepDep_FV;
  ///\brief [GEANT4]: The total 'early' lepton energy deposited within the stops
  /// veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. Unlike other branches, this does *not* perform descendent roll up.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double LepDep_veto;
  ///\brief [GEANT4]: The total 'early' energy deposited within the stops
  /// fiducial volume by descendents of the primary lepton.
  ///
  /// N.B. This branch is most useful for determining energy deposited by
  /// primary muon descendents, which will likely be michel electrons.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double LepDepDescendent_FV;
  ///\brief [GEANT4]: The total 'early' lepton energy deposited within the stops
  /// veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch is most useful for determining energy deposited by
  /// primary muon descendents, which will likely be michel electrons.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double LepDepDescendent_veto;
  ///\brief [GEANT4]: The total 'early' proton energy deposited within the stops
  /// fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double ProtonDep_FV;
  ///\brief [GEANT4]: The total 'early' proton energy deposited within the stops
  /// veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double ProtonDep_veto;
  ///\brief [GEANT4]: The total 'early' neutron energy deposited within the
  /// stops fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double NeutronDep_FV;
  ///\brief [GEANT4]: The charge-weighted average time of all neutron deposites
  /// within the stops fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double NeutronDep_ChrgWAvgTime_FV;
  ///\brief [GEANT4]: The total 'early' neutron energy deposited within the
  /// stops veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double NeutronDep_veto;
  ///\brief [GEANT4]: The charge-weighted average time of all neutron deposites
  /// within the stops veto volume, but within the active LAr volume of the
  /// stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double NeutronDep_ChrgWAvgTime_veto;
  ///\brief [GEANT4]: The total 'early' charged pion energy deposited within the
  /// stops fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double PiCDep_FV;
  ///\brief [GEANT4]: The total 'early' charged pion energy deposited within the
  /// stops veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double PiCDep_veto;
  ///\brief [GEANT4]: The total 'early' neutral pion energy deposited within the
  /// stops fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double Pi0Dep_FV;
  ///\brief [GEANT4]: The total 'early' neutral pion energy deposited within the
  /// stops veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double Pi0Dep_veto;
  ///\brief [GEANT4]: The total 'early' 'other' particle energy deposited within
  /// the stops fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double OtherDep_FV;
  ///\brief [GEANT4]: The total 'early' 'other' particle energy deposited within
  /// the stops veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double OtherDep_veto;

  ///\brief [GEANT4]: The total 'early' non-GENIE-simulated-lepton particle
  /// energy deposited within thestops veto volume, but within the active LAr
  /// volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double TotalNonlep_Dep_FV;
  ///\brief [GEANT4]: The total 'early' non-GENIE-simulated-lepton particle
  /// energy deposited within thestops veto volume, but within the active LAr
  /// volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *before* the time separator, if none was used these
  /// contain the energy integrated over all simulation time
  double TotalNonlep_Dep_veto;

  ///\brief [GEANT4]: The total 'late' lepton energy deposited within the stops
  /// fiducial volume.
  ///
  /// N.B. Unlike other branches, this does *not* perform descendent roll up.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double LepDep_timesep_FV;
  ///\brief [GEANT4]: The total 'late' lepton energy deposited within the stops
  /// veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. Unlike other branches, this does *not* perform descendent roll up.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double LepDep_timesep_veto;
  ///\brief [GEANT4]: The total 'late' energy deposited within the stops
  /// fiducial volume by descendents of the primary lepton.
  ///
  /// N.B. This branch is most useful for determining energy deposited by
  /// primary muon descendents, which will likely be michel electrons.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double LepDepDescendent_timesep_FV;
  ///\brief [GEANT4]: The total 'late' lepton energy deposited within the stops
  /// veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch is most useful for determining energy deposited by
  /// primary muon descendents, which will likely be michel electrons.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double LepDepDescendent_timesep_veto;
  ///\brief [GEANT4]: The total 'late' proton energy deposited within the stops
  /// fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double ProtonDep_timesep_FV;
  ///\brief [GEANT4]: The total 'late' proton energy deposited within the stops
  /// veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double ProtonDep_timesep_veto;
  ///\brief [GEANT4]: The total 'late' neutron energy deposited within the
  /// stops fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double NeutronDep_timesep_FV;
  ///\brief [GEANT4]: The total 'late' neutron energy deposited within the
  /// stops veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double NeutronDep_timesep_veto;
  ///\brief [GEANT4]: The total 'late' charged pion energy deposited within the
  /// stops fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double PiCDep_timesep_FV;
  ///\brief [GEANT4]: The total 'late' charged pion energy deposited within the
  /// stops veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double PiCDep_timesep_veto;
  ///\brief [GEANT4]: The total 'late' neutral pion energy deposited within the
  /// stops fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double Pi0Dep_timesep_FV;
  ///\brief [GEANT4]: The total 'late' neutral pion energy deposited within the
  /// stops veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double Pi0Dep_timesep_veto;
  ///\brief [GEANT4]: The total 'late' 'other' particle energy deposited within
  /// the stops fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double OtherDep_timesep_FV;
  ///\brief [GEANT4]: The total 'late' 'other' particle energy deposited within
  /// the stops veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double OtherDep_timesep_veto;

  ///\brief [GEANT4]: The total 'late' non-GENIE-simulated-lepton particle
  /// energy deposited within thestops veto volume, but within the active LAr
  /// volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double TotalNonlep_Dep_timesep_FV;
  ///\brief [GEANT4]: The total 'late' non-GENIE-simulated-lepton particle
  /// energy deposited within thestops veto volume, but within the active LAr
  /// volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  ///
  /// N.B. If this was run with a deposit time separator, these branches contain
  /// the energy deposited *after* the time separator, if none was used these
  /// will not be filled.
  double TotalNonlep_Dep_timesep_veto;

  ///\brief [GEANT4]: Whether the primary lepton left the active stop volume.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExit;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop volume
  /// with more KE than a runtime threshold (default = 50 MeV);
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExit_AboveThresh;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop via the +Z
  /// face.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExitBack;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop via the -Z
  /// face.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExitFront;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop via the -Y
  /// face.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExitYLow;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop via the +Y
  /// face.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExitYHigh;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop via the -X
  /// face.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExitXLow;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop via the +X
  /// face.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExitXHigh;
  ///\brief [GEANT4]: The exit topology of the primary lepton.
  ///
  /// *  0: Did not exit
  /// *  1: Exit Back
  /// *  2: Exit Front
  /// *  3: Exit Y Low
  /// *  4: Exit Y High
  /// *  5: Exit X Low
  /// *  6: Exit X High
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  int LepExitTopology;

  /// [GEANT4]: The exit KE of the primary lepton.
  double LepExitKE;
  /// [GEANT4]: The exit 3-position of the primary lepton.
  double LepExitingPos[3];
  /// [GEANT4]: The exit 3-momentum of the primary lepton.
  double LepExitingMom[3];

  ///\brief [EVENT SUMMARY]: Whether interaction involved a (anti-) muon
  /// neutrino
  bool IsNumu;
  ///\brief [EVENT SUMMARY]: Whether interaction an anti-neutrino
  bool IsAntinu;
  ///\brief [EVENT SUMMARY]: Whether interaction was charged current
  bool IsCC;
  ///\brief [EVENT SUMMARY]: Whether the GENIE simulation produced no final
  /// state pions.
  bool Is0Pi;
  ///\brief [EVENT SUMMARY]: Whether the GENIE simulation produced one final
  /// state charged pion.
  bool Is1PiC;
  ///\brief [EVENT SUMMARY]: Whether the GENIE simulation produced one final
  /// state neutral pion.
  bool Is1Pi0;
  ///\brief [EVENT SUMMARY]: Whether the GENIE simulation produced one final
  /// state pion.
  bool Is1Pi;
  ///\brief [EVENT SUMMARY]: Whether the GENIE simulation produced multiple
  /// final state pions.
  bool IsNPi;
  ///\brief [EVENT SUMMARY]: Whether the GENIE simulation produced other final
  /// state particles.
  ///
  /// N.B. This is often due to gamma or kaon emission.
  bool IsOther;
  ///\brief [EVENT SUMMARY]: The summarised event topology
  ///
  /// Negative numbers indicate NC interactions.
  /// * 1 : 0Pi
  /// * 2 : 1PiC
  /// * 3 : 1Pi0
  /// * 4 : NPi
  /// * 5 : other
  int Topology;
  ///\brief [EVENT SUMMARY]: Whether the hadronic shower is contained within the
  /// stop fiducial volume.
  ///
  /// N.B. This checks whether the total veto-region deposit is greater than the
  /// threshold passed by command line (or 10 MeV by default.). This can be
  /// fully recalculated given a different threshold by summing over the
  /// XXXXDep_veto branches.
  bool HadrShowerContainedInFV;
  ///\brief [EVENT SUMMARY]: Whether the primary lepton deposits are contained
  /// within the stop fiducial volume.
  ///
  /// N.B. This is useful for checking whether electron neutrino events had
  /// contain EM showers, it is less useful for muon neutrino interactions.
  ///
  /// N.B. This checks whether the total veto-region deposit is greater than the
  /// threshold passed by command line (or 10 MeV by default.). This can be
  /// fully recalculated given a different threshold by summing over the
  /// LepDep_veto branch.
  bool PrimaryLeptonContainedInFV;

  void Reset();

  void Copy(DepositsSummary const &other);

  static DepositsSummary *MakeTreeWriter(TTree *OutputTree, double timesep_us = 0xdeadbeef,
                              bool IsLite = false);
  void SetBranchAddresses();
};

#endif
