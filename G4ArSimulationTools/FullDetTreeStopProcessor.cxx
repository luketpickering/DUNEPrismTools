#include "DetectorStop.hxx"
#include "Utils.hxx"

#include "FullDetTreeReader.h"
#include "G4ArReader.h"

#include "TH3D.h"
#include "TRandom.h"
#include "TTree.h"

#include <utility>

struct DetBox {
  double XOffset;
  double XWidth_fv;
  double YWidth_fv;
  double ZWidth_fv;

  double XWidth_det;
  double YWidth_det;
  double ZWidth_det;

  double X_Range_fv[2];
  double Y_Range_fv[2];
  double Z_Range_fv[2];

  size_t X_fv[2];
  size_t X_veto_left[2];
  size_t X_veto_right[2];

  double POTExposure;
};

/// Energy deposit and GENIE passthrough output tree
struct EDep {
  ///\brief The detector stop number used, refer to input xml for stop offsets.
  ///
  /// N.B. When overlapping stops are defined the event is randomly placed
  /// within one of the overlapping stops at the interaction position. The
  /// choice is weighted by the POTExposure branch in the input run plan
  /// xml.
  int stop;

  /// [GENIE P/T]:  The vertex 3-position in cm
  double vtx[3];
  ///\brief [GENIE P/T]:  The X position of the vertex relative to the centre of
  ///a
  /// stop in cm.
  double vtxInDetX;
  /// [GENIE P/T]:  The X offset of stop, stop in cm.
  double XOffset;
  /// [GENIE P/T]:  The GENIE interaction code.
  TObjString *EventCode;
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
  ///\brief [GENIE P/T]:  The number of final state other particles in the
  ///event.
  ///
  /// N.B. These do not include GENIE bindinos or nuclear PDG codes.
  /// By eye, these are most often Kaons or Lambdas.
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
  ///\brief [GENIE P/T]: The total energy of all non-primary leptons at the
  /// end of the GENIE simulation.
  double Total_ENonPrimaryLep_True;

  ///\brief [GEANT4]: The total lepton energy deposited within the stops
  /// fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation. i.e. This quantity will likely contain deposits
  /// from Michel electrons from stopped primary muon decays.
  double LepDep_FV;
  ///\brief [GEANT4]: The total lepton energy deposited within the stops
  /// veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation. i.e. This quantity will likely contain deposits
  /// from Michel electrons from stopped primary muon decays.
  double LepDep_veto;
  ///\brief [GEANT4]: The total proton energy deposited within the stops
  /// fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double ProtonDep_FV;
  ///\brief [GEANT4]: The total proton energy deposited within the stops
  /// veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double ProtonDep_veto;
  ///\brief [GEANT4]: The total neutron energy deposited within the stops
  /// fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double NeutronDep_FV;
  ///\brief [GEANT4]: The charge-weighted average time of all neutron deposites
  /// within the stops fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double NeutronDep_ChrgWAvgTime_FV;
  ///\brief [GEANT4]: The total neutron energy deposited within the stops
  /// veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double NeutronDep_veto;
  ///\brief [GEANT4]: The charge-weighted average time of all neutron deposites
  /// within the stops veto volume, but within the active LAr volume of the
  /// stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double NeutronDep_ChrgWAvgTime_veto;
  ///\brief [GEANT4]: The total charged pion energy deposited within the stops
  /// fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double PiCDep_FV;
  ///\brief [GEANT4]: The total charged pion energy deposited within the stops
  /// veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double PiCDep_veto;
  ///\brief [GEANT4]: The total neutral pion energy deposited within the stops
  /// fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double Pi0Dep_FV;
  ///\brief [GEANT4]: The total neutral pion energy deposited within the stops
  /// veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double Pi0Dep_veto;
  ///\brief [GEANT4]: The total 'other' particle energy deposited within the
  ///  stops fiducial volume.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double OtherDep_FV;
  ///\brief [GEANT4]: The total 'other' particle energy deposited within the
  ///  stops veto volume, but within the active LAr volume of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double OtherDep_veto;

  ///\brief [GEANT4]: The total non-GENIE-simulated-lepton particle energy
  /// deposited within thestops veto volume, but within the active LAr volume
  /// of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double TotalNonlep_Dep_FV;
  ///\brief [GEANT4]: The total non-GENIE-simulated-lepton particle energy
  /// deposited within thestops veto volume, but within the active LAr volume
  /// of the stop.
  ///
  /// N.B. This branch rolls up all deposits by all descendent particles in
  /// the GEANT4 simulation.
  double TotalNonlep_Dep_veto;

  ///\brief [GEANT4]: Whether the primary lepton left the active stop volume.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExit;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop via the +Z
  ///face.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExitBack;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop via the -Z
  ///face.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExitFront;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop via the -Y
  ///face.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExitYLow;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop via the +Y
  ///face.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExitYHigh;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop via the -X
  ///face.
  ///
  /// N.B. This will track a primary electron, but that should shower very
  /// quickly. This branch is nominally designed for primary muons.
  bool LepExitXLow;
  ///\brief [GEANT4]: Whether the primary lepton left the active stop via the +X
  ///face.
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

  /// [GEANT4]: The exit 3-position of the primary lepton.
  double LepExitingPos[3];
  /// [GEANT4]: The exit 3-momentum of the primary lepton.
  double LepExitingMom[3];

  ///\brief [EVENT SUMMARY]: Whether interaction involved a (anti-) muon
  ///neutrino
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
  /// fully recalculated given a different threshold by summing over the the
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
  /// fully recalculated given a different threshold by summing over the the
  /// LepDep_veto branch.
  bool PrimaryLeptonContainedInFV;
};

std::vector<DetectorStop> DetectorStops;
std::vector<DetBox> Detectors;

std::string inpDir, outputFile;
std::string runPlanCfg, runPlanName = "";
double VetoThreshold = 0.001;

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "\t-i <input dir>                : Input directory containing "
         "condensed arbox.py output.\n"
         "\t-r <RunPlan.XML>              : An XML file specifying a run "
         "plan  \n"
         "\t                                to build fluxes for. See     "
         "      \n"
         "\t                                documentation for XML "
         "structure.   \n"
         "\t-o <output.root>              : Output file name.\n"
         "\t-v <veto threshold>           : Threshold energy deposit in "
         "veto region to pass \n"
         "\t                                selection {default = 10 MeV}.\n"
      << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-?") {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      inpDir = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      outputFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-r") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");

      runPlanCfg = params[0];

      if (params.size() > 1) {
        runPlanName = params[1];
      }
    } else if (std::string(argv[opt]) == "-v") {
      VetoThreshold = str2T<double>(argv[++opt]);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

enum OOFVInclusion {
  kIncludeFVYZ = 1,
  kIncludeOOFVYZ = 2,
  kIncludeWholeDet = 3
};
double Jaccumulate(double ***arr, size_t ind1, size_t ind2, double init,
                   OOFVInclusion fvincl) {
  double sum = init;
  for (size_t i = ind1; i < ind2; ++i) {
    if ((arr[i][1][1] != 0) && (!std::isnormal(arr[i][1][1]))) {
      std::cout << "[" << i << "][1][1]: bin non-normal = " << arr[i][1][1]
                << std::endl;
      throw;
    }

    if (fvincl & 1) {
      sum += arr[i][1][1];
    }

    if (fvincl & 2) {
      sum += arr[i][0][0];
      sum += arr[i][0][1];
      sum += arr[i][0][2];

      sum += arr[i][1][0];
      sum += arr[i][1][2];

      sum += arr[i][2][0];
      sum += arr[i][2][1];
      sum += arr[i][2][2];

      if ((arr[i][0][0] != 0) && (!std::isnormal(arr[i][0][0]))) {
        std::cout << "[" << i << "][0][0]: bin non-normal = " << arr[i][0][0]
                  << std::endl;
        throw;
      }
      if ((arr[i][0][1] != 0) && (!std::isnormal(arr[i][0][1]))) {
        std::cout << "[" << i << "][0][1]: bin non-normal = " << arr[i][0][1]
                  << std::endl;
        throw;
      }
      if ((arr[i][0][2] != 0) && (!std::isnormal(arr[i][0][2]))) {
        std::cout << "[" << i << "][0][2]: bin non-normal = " << arr[i][0][2]
                  << std::endl;
        throw;
      }

      if ((arr[i][1][0] != 0) && (!std::isnormal(arr[i][1][0]))) {
        std::cout << "[" << i << "][1][0]: bin non-normal = " << arr[i][1][0]
                  << std::endl;
        throw;
      }
      if ((arr[i][1][2] != 0) && (!std::isnormal(arr[i][1][2]))) {
        std::cout << "[" << i << "][1][2]: bin non-normal = " << arr[i][1][2]
                  << std::endl;
        throw;
      }

      if ((arr[i][2][0] != 0) && (!std::isnormal(arr[i][2][0]))) {
        std::cout << "[" << i << "][2][0]: bin non-normal = " << arr[i][2][0]
                  << std::endl;
        throw;
      }
      if ((arr[i][2][1] != 0) && (!std::isnormal(arr[i][2][1]))) {
        std::cout << "[" << i << "][2][1]: bin non-normal = " << arr[i][2][1]
                  << std::endl;
        throw;
      }
      if ((arr[i][2][2] != 0) && (!std::isnormal(arr[i][2][2]))) {
        std::cout << "[" << i << "][2][2]: bin non-normal = " << arr[i][2][2]
                  << std::endl;
        throw;
      }
    }
  }
  return sum;
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  TChain *config = new TChain("configTree");
  config->Add(inpDir.c_str());

  DetectorAndFVDimensions detdims;
  Int_t NMaxTrackSteps;

  config->SetBranchAddress("NXSteps", &detdims.NXSteps);
  config->SetBranchAddress("DetMin", &detdims.DetMin);
  config->SetBranchAddress("DetMax", &detdims.DetMax);
  config->SetBranchAddress("FVGap", &detdims.FVGap);
  config->SetBranchAddress("NMaxTrackSteps", &NMaxTrackSteps);

  config->GetEntry(0);

  TH3D *DetMap = detdims.BuildDetectorMap();

  FullDetTreeReader *rdr = new FullDetTreeReader(
      "fulldetTree", inpDir, detdims.NXSteps, NMaxTrackSteps);

  if (!rdr->GetEntries()) {
    std::cout << "No valid input files found." << std::endl;
    return 1;
  }
  if (!runPlanCfg.length()) {
    std::cout << "[ERROR]: Found no run plan configuration file." << std::endl;
  }

  TFile *outfile = CheckOpenFile(outputFile, "RECREATE");

  // Read in det stops
  DetectorStops = ReadDetectorStopConfig(runPlanCfg, runPlanName);

  size_t NDets = DetectorStops.size();
  Detectors.reserve(NDets);

  EDep OutputEDep;
  TTree *OutputTree = new TTree("EDeps", "");

  // translate to detboxes
  for (size_t d_it = 0; d_it < NDets; ++d_it) {
    DetBox db;

    db.XOffset = -DetectorStops[d_it].LateralOffset * 100.0;
    db.XWidth_fv = DetectorStops[d_it].DetectorFiducialWidth * 100.0;
    db.YWidth_fv = DetectorStops[d_it].DetectorFiducialHeight * 100.0;
    db.ZWidth_fv = DetectorStops[d_it].DetectorFiducialDepth * 100.0;
    db.XWidth_det = db.XWidth_fv + 2 * detdims.FVGap[0];
    db.YWidth_det = db.YWidth_fv + 2 * detdims.FVGap[1];
    db.ZWidth_det = db.ZWidth_fv + 2 * detdims.FVGap[2];
    db.POTExposure = DetectorStops[d_it].POTExposure;

    double Detlow = db.XOffset - db.XWidth_det / 2.0;
    double DetHigh = db.XOffset + db.XWidth_det / 2.0;

    db.X_Range_fv[0] = db.XOffset - db.XWidth_fv / 2.0;
    db.X_Range_fv[1] = db.XOffset + db.XWidth_fv / 2.0;

    db.Y_Range_fv[0] = -db.YWidth_fv / 2.0;
    db.Y_Range_fv[1] = db.YWidth_fv / 2.0;

    db.Z_Range_fv[0] = -db.ZWidth_fv / 2.0;
    db.Z_Range_fv[1] = db.ZWidth_fv / 2.0;

    std::cout << "[INFO]: Det: [" << Detlow << ", " << DetHigh << "], FV: ["
              << db.X_Range_fv[0] << ", " << db.X_Range_fv[1] << "]."
              << std::endl;

    db.X_fv[0] = DetMap->GetXaxis()->FindFixBin(db.X_Range_fv[0]) - 1;
    db.X_fv[1] = DetMap->GetXaxis()->FindFixBin(db.X_Range_fv[1]) - 1;

    db.X_veto_left[0] = DetMap->GetXaxis()->FindFixBin(Detlow) - 1;
    db.X_veto_left[1] =
        DetMap->GetXaxis()->FindFixBin(Detlow + detdims.FVGap[0]) - 1;

    db.X_veto_right[0] =
        DetMap->GetXaxis()->FindFixBin(DetHigh - detdims.FVGap[0]) - 1;
    db.X_veto_right[1] = DetMap->GetXaxis()->FindFixBin(DetHigh) - 1;

    Detectors.push_back(db);
  }

  OutputTree->Branch("stop", &OutputEDep.stop, "stop/I");

  OutputTree->Branch("vtx", &OutputEDep.vtx, "vtx[3]/D");
  OutputTree->Branch("vtxInDetX", &OutputEDep.vtxInDetX, "vtxInDetX/D");
  OutputTree->Branch("XOffset", &OutputEDep.XOffset, "XOffset/D");
  OutputEDep.EventCode = nullptr;
  OutputTree->Branch("EventCode", &OutputEDep.EventCode);

  OutputTree->Branch("nu_4mom", &OutputEDep.nu_4mom, "nu_4mom[4]/D");
  OutputTree->Branch("y_True", &OutputEDep.y_True, "y_True/D");
  OutputTree->Branch("W_Rest", &OutputEDep.W_Rest, "W_Rest/D");
  OutputTree->Branch("Q2_True", &OutputEDep.Q2_True, "Q2_True/D");
  OutputTree->Branch("FourMomTransfer_True", &OutputEDep.FourMomTransfer_True,
                     "FourMomTransfer_True[4]/D");

  OutputTree->Branch("nu_PDG", &OutputEDep.nu_PDG, "nu_PDG/I");
  OutputTree->Branch("PrimaryLepPDG", &OutputEDep.PrimaryLepPDG,
                     "PrimaryLepPDG/I");
  OutputTree->Branch("PrimaryLep_4mom", &OutputEDep.PrimaryLep_4mom,
                     "PrimaryLep_4mom[4]/D");

  OutputTree->Branch("NLep", &OutputEDep.NLep, "NLep/I");
  OutputTree->Branch("NPi0", &OutputEDep.NPi0, "NPi0/I");
  OutputTree->Branch("NPiC", &OutputEDep.NPiC, "NPiC/I");
  OutputTree->Branch("NProton", &OutputEDep.NProton, "NProton/I");
  OutputTree->Branch("NNeutron", &OutputEDep.NNeutron, "NNeutron/I");
  OutputTree->Branch("NGamma", &OutputEDep.NGamma, "NGamma/I");
  OutputTree->Branch("NOther", &OutputEDep.NOther, "NOther/I");

  OutputTree->Branch("EKinPi0_True", &OutputEDep.EKinPi0_True,
                     "EKinPi0_True/D");
  OutputTree->Branch("EMassPi0_True", &OutputEDep.EMassPi0_True,
                     "EMassPi0_True/D");
  OutputTree->Branch("EKinPiC_True", &OutputEDep.EKinPiC_True,
                     "EKinPiC_True/D");
  OutputTree->Branch("EMassPiC_True", &OutputEDep.EMassPiC_True,
                     "EMassPiC_True/D");
  OutputTree->Branch("EKinProton_True", &OutputEDep.EKinProton_True,
                     "EKinProton_True/D");
  OutputTree->Branch("EMassProton_True", &OutputEDep.EMassProton_True,
                     "EMassProton_True/D");
  OutputTree->Branch("EKinNeutron_True", &OutputEDep.EKinNeutron_True,
                     "EKinNeutron_True/D");
  OutputTree->Branch("EMassNeutron_True", &OutputEDep.EMassNeutron_True,
                     "EMassNeutron_True/D");
  OutputTree->Branch("EGamma_True", &OutputEDep.EGamma_True, "EGamma_True/D");
  OutputTree->Branch("EOther_True", &OutputEDep.EOther_True, "EOther_True/D");

  OutputTree->Branch("Total_ENonPrimaryLep_True",
                     &OutputEDep.Total_ENonPrimaryLep_True,
                     "Total_ENonPrimaryLep_True/D");

  OutputTree->Branch("LepDep_FV", &OutputEDep.LepDep_FV, "LepDep_FV/D");
  OutputTree->Branch("LepDep_veto", &OutputEDep.LepDep_veto, "LepDep_veto/D");
  OutputTree->Branch("ProtonDep_FV", &OutputEDep.ProtonDep_FV,
                     "ProtonDep_FV/D");
  OutputTree->Branch("ProtonDep_veto", &OutputEDep.ProtonDep_veto,
                     "ProtonDep_veto/D");
  OutputTree->Branch("NeutronDep_FV", &OutputEDep.NeutronDep_FV,
                     "NeutronDep_FV/D");
  OutputTree->Branch("NeutronDep_ChrgWAvgTime_FV",
                     &OutputEDep.NeutronDep_ChrgWAvgTime_FV,
                     "NeutronDep_ChrgWAvgTime_FV/D");
  OutputTree->Branch("NeutronDep_veto", &OutputEDep.NeutronDep_veto,
                     "NeutronDep_veto/D");
  OutputTree->Branch("NeutronDep_ChrgWAvgTime_veto",
                     &OutputEDep.NeutronDep_ChrgWAvgTime_veto,
                     "NeutronDep_ChrgWAvgTime_veto/D");
  OutputTree->Branch("PiCDep_FV", &OutputEDep.PiCDep_FV, "PiCDep_FV/D");
  OutputTree->Branch("PiCDep_veto", &OutputEDep.PiCDep_veto, "PiCDep_veto/D");
  OutputTree->Branch("Pi0Dep_FV", &OutputEDep.Pi0Dep_FV, "Pi0Dep_FV/D");
  OutputTree->Branch("Pi0Dep_veto", &OutputEDep.Pi0Dep_veto, "Pi0Dep_veto/D");
  OutputTree->Branch("OtherDep_FV", &OutputEDep.OtherDep_FV, "OtherDep_FV/D");
  OutputTree->Branch("OtherDep_veto", &OutputEDep.OtherDep_veto,
                     "OtherDep_veto/D");

  OutputTree->Branch("TotalNonlep_Dep_FV", &OutputEDep.TotalNonlep_Dep_FV,
                     "TotalNonlep_Dep_FV/D");
  OutputTree->Branch("TotalNonlep_Dep_veto", &OutputEDep.TotalNonlep_Dep_veto,
                     "TotalNonlep_Dep_veto/D");

  OutputTree->Branch("LepExit", &OutputEDep.LepExit, "LepExit/O");

  OutputTree->Branch("LepExitBack", &OutputEDep.LepExitBack, "LepExitBack/O");
  OutputTree->Branch("LepExitFront", &OutputEDep.LepExitFront,
                     "LepExitFront/O");
  OutputTree->Branch("LepExitYLow", &OutputEDep.LepExitYLow, "LepExitYLow/O");
  OutputTree->Branch("LepExitYHigh", &OutputEDep.LepExitYHigh,
                     "LepExitYHigh/O");
  OutputTree->Branch("LepExitXLow", &OutputEDep.LepExitXLow, "LepExitXLow/O");
  OutputTree->Branch("LepExitXHigh", &OutputEDep.LepExitXHigh,
                     "LepExitXHigh/O");

  OutputTree->Branch("LepExitTopology", &OutputEDep.LepExitTopology,
                     "LepExitTopology/I");

  OutputTree->Branch("LepExitingPos", &OutputEDep.LepExitingPos,
                     "LepExitingPos[3]/D");
  OutputTree->Branch("LepExitingMom", &OutputEDep.LepExitingMom,
                     "LepExitingMom[3]/D");

  OutputTree->Branch("IsNumu", &OutputEDep.IsNumu, "IsNumu/O");
  OutputTree->Branch("IsAntinu", &OutputEDep.IsAntinu, "IsAntinu/O");
  OutputTree->Branch("IsCC", &OutputEDep.IsCC, "IsCC/O");
  OutputTree->Branch("Is0Pi", &OutputEDep.Is0Pi, "Is0Pi/O");
  OutputTree->Branch("Is1PiC", &OutputEDep.Is1PiC, "Is1PiC/O");
  OutputTree->Branch("Is1Pi0", &OutputEDep.Is1Pi0, "Is1Pi0/O");
  OutputTree->Branch("Is1Pi", &OutputEDep.Is1Pi, "Is1Pi/O");
  OutputTree->Branch("IsNPi", &OutputEDep.IsNPi, "IsNPi/O");
  OutputTree->Branch("IsOther", &OutputEDep.IsOther, "IsOther/O");
  OutputTree->Branch("Topology", &OutputEDep.Topology, "Topology/I");
  OutputTree->Branch("HadrShowerContainedInFV",
                     &OutputEDep.HadrShowerContainedInFV,
                     "HadrShowerContainedInFV/O");
  OutputTree->Branch("PrimaryLeptonContainedInFV",
                     &OutputEDep.PrimaryLeptonContainedInFV,
                     "PrimaryLeptonContainedInFV/O");

  std::cout << "[INFO]: Reading " << rdr->GetEntries() << " input entries."
            << std::endl;

  size_t loud_every = rdr->GetEntries() / 10;
  TRandom *rand = new TRandom();

  std::vector<int> overlapping_stops;
  for (size_t e_it = 0; e_it < rdr->GetEntries(); ++e_it) {
    rdr->GetEntry(e_it);

    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << rdr->VertexPosition[0] << ", " << rdr->VertexPosition[1]
                << ", " << rdr->VertexPosition[2]
                << "}, Enu: " << rdr->nu_4mom[3] << " )" << std::flush;
    }

    // DetBox + Reader -> EDep -> Fill out tree
    //
    //
    // Jake: Adding in a vector of overlapping stops

    overlapping_stops.clear();

    int stop = -1;
    for (size_t d_it = 0; d_it < NDets; ++d_it) {
      DetBox &db = Detectors[d_it];

      if ((rdr->VertexPosition[0] < db.X_Range_fv[0]) ||
          (rdr->VertexPosition[0] > db.X_Range_fv[1]) ||
          (rdr->VertexPosition[1] < db.Y_Range_fv[0]) ||
          (rdr->VertexPosition[1] > db.Y_Range_fv[1]) ||
          (rdr->VertexPosition[2] < db.Z_Range_fv[0]) ||
          (rdr->VertexPosition[2] > db.Z_Range_fv[1])) {
        continue;
      }
      stop = d_it;
      overlapping_stops.push_back(stop);
    }

    if (stop == -1) {
      continue;
    }

    if (overlapping_stops.size() > 1) {  // Choose stop according to exposure

      double total_POT = 0.;

      for (size_t s_it = 0; s_it < overlapping_stops.size(); ++s_it) {
        total_POT += Detectors.at(s_it).POTExposure;
      }

      double rand_val = rand->Uniform() * total_POT;
      double running_POT = 0;

      stop = overlapping_stops.back();
      for (size_t s_it = 0; s_it < overlapping_stops.size(); ++s_it) {
        running_POT += Detectors.at(s_it).POTExposure;
        if (rand_val < running_POT) {
          stop = overlapping_stops[s_it];
          break;
        }
      }
    }

    DetBox &stopBox = Detectors[stop];
    OutputEDep.stop = stop;
    (*OutputEDep.EventCode) = (*rdr->EventCode);
    std::copy_n(rdr->VertexPosition, 3, OutputEDep.vtx);

    OutputEDep.vtxInDetX = OutputEDep.vtx[0] - stopBox.XOffset;
    OutputEDep.XOffset = stopBox.XOffset;

    if (fabs(OutputEDep.vtxInDetX) > stopBox.XWidth_det / 2.0) {
      std::cout << OutputEDep.vtxInDetX << std::endl;
      throw;
    }

    std::copy_n(rdr->nu_4mom, 4, OutputEDep.nu_4mom);
    OutputEDep.nu_PDG = rdr->nu_PDG;

    OutputEDep.y_True = rdr->y_True;
    OutputEDep.Q2_True = rdr->Q2_True;
    std::copy_n(rdr->FourMomTransfer_True, 4, OutputEDep.FourMomTransfer_True);
    OutputEDep.W_Rest = rdr->W_Rest;

    OutputEDep.PrimaryLepPDG = rdr->PrimaryLepPDG;
    std::copy_n(rdr->PrimaryLep_4mom, 4, OutputEDep.PrimaryLep_4mom);

    OutputEDep.NLep = rdr->NLep;
    OutputEDep.NPi0 = rdr->NPi0;
    OutputEDep.NPiC = rdr->NPiC;
    OutputEDep.NProton = rdr->NProton;
    OutputEDep.NNeutron = rdr->NNeutron;
    OutputEDep.NGamma = rdr->NGamma;
    OutputEDep.NOther = rdr->NOther;

    OutputEDep.EKinPi0_True = rdr->EKinPi0_True;
    OutputEDep.EMassPi0_True = rdr->EMassPi0_True;
    OutputEDep.EKinPiC_True = rdr->EKinPiC_True;
    OutputEDep.EMassPiC_True = rdr->EMassPiC_True;
    OutputEDep.EKinProton_True = rdr->EKinProton_True;
    OutputEDep.EMassProton_True = rdr->EMassProton_True;
    OutputEDep.EKinNeutron_True = rdr->EKinNeutron_True;
    OutputEDep.EMassNeutron_True = rdr->EMassNeutron_True;
    OutputEDep.EGamma_True = rdr->EGamma_True;
    OutputEDep.Total_ENonPrimaryLep_True = rdr->ENonPrimaryLep_True;

    double DetXLow = stopBox.XOffset - stopBox.XWidth_det / 2.0;
    double DetXHigh = stopBox.XOffset + stopBox.XWidth_det / 2.0;
    double DetYLow = 0. - stopBox.YWidth_det / 2.0;
    double DetYHigh = stopBox.YWidth_det / 2.0;
    double DetZLow = 0. - stopBox.ZWidth_det / 2.0;
    double DetZHigh = stopBox.ZWidth_det / 2.0;

    // Checking if lepton exits -- muon only
    OutputEDep.LepExitBack = false;
    OutputEDep.LepExitFront = false;
    OutputEDep.LepExitYLow = false;
    OutputEDep.LepExitYHigh = false;
    OutputEDep.LepExitXLow = false;
    OutputEDep.LepExitXHigh = false;
    OutputEDep.LepExit = false;

    if (abs(rdr->PrimaryLepPDG) == 13) {
      for (int i = 0; i < rdr->NMuonTrackSteps; ++i) {
        if (OutputEDep.LepExit) {
          break;
        }

        double posX = rdr->MuonTrackPos[i][0];
        double posY = rdr->MuonTrackPos[i][1];
        double posZ = rdr->MuonTrackPos[i][2];
        if ((posX == 0xdeadbeef) && (posY == 0xdeadbeef) &&
            (posZ == 0xdeadbeef)) {
          std::cout << "[INFO]: Step " << i << ", lep pos deadbeef."
                    << std::endl;
          break;
        }

        double dX = 0., dY = 0., dZ = 0.;
        bool exitX = false, exitY = false, exitZ = false;

        if (posX <= DetXLow || posX >= DetXHigh) {
          dX = fabs(DetXLow - posX);
          exitX = true;
        } else {
          exitX = false;
        }

        if (posY <= DetYLow || posY >= DetYHigh) {
          dY = fabs(DetYLow - posY);
          exitY = true;
        } else {
          exitY = false;
        }

        if (posZ <= DetZLow || posZ >= DetZHigh) {
          dZ = fabs(DetZLow - posZ);
          exitZ = true;
        } else {
          exitZ = false;
        }

        if (exitX || exitY || exitZ) {
          if (dX > dY && dX > dZ) {
            OutputEDep.LepExitBack = false;
            OutputEDep.LepExitFront = false;
            OutputEDep.LepExitYLow = false;
            OutputEDep.LepExitYHigh = false;

            if (posX <= DetXLow) {
              OutputEDep.LepExitXLow = true;
              OutputEDep.LepExitXHigh = false;
            } else if (posX >= DetXHigh) {
              OutputEDep.LepExitXLow = false;
              OutputEDep.LepExitXHigh = true;
            } else {
              std::cout << "ERROR X" << std::endl;
              throw;
            }
          } else if (dY > dX && dY > dZ) {
            OutputEDep.LepExitBack = false;
            OutputEDep.LepExitFront = false;
            OutputEDep.LepExitXLow = false;
            OutputEDep.LepExitXHigh = false;

            if (posY <= DetYLow) {
              OutputEDep.LepExitYLow = true;
              OutputEDep.LepExitYHigh = false;
            } else if (posY >= DetYHigh) {
              OutputEDep.LepExitYLow = false;
              OutputEDep.LepExitYHigh = true;
            } else {
              std::cout << "ERROR Y" << std::endl;
              throw;
            }
          } else if (dZ > dX && dZ > dY) {
            OutputEDep.LepExitXLow = false;
            OutputEDep.LepExitXHigh = false;
            OutputEDep.LepExitYLow = false;
            OutputEDep.LepExitYHigh = false;

            if (posZ <= DetZLow) {
              OutputEDep.LepExitBack = false;
              OutputEDep.LepExitFront = true;
            } else if (posZ >= DetZHigh) {
              OutputEDep.LepExitFront = false;
              OutputEDep.LepExitBack = true;
            } else {
              std::cout << "ERROR Z" << std::endl;
              throw;
            }
          }
          OutputEDep.LepExit = true;
          std::copy_n(rdr->MuonTrackPos[i], 3, OutputEDep.LepExitingPos);
          std::copy_n(rdr->MuonTrackMom[i], 3, OutputEDep.LepExitingMom);

#ifdef DEBUG
          std::cout << "[INFO]: Muon stopped tracking at step [" << i
                    << "] at {" << rdr->MuonTrackPos[i][0] << ", "
                    << rdr->MuonTrackPos[i][1] << ", "
                    << rdr->MuonTrackPos[i][2] << " } with 3-mom {"
                    << rdr->MuonTrackMom[i][0] << ", "
                    << rdr->MuonTrackMom[i][1] << ", "
                    << rdr->MuonTrackMom[i][2] << " }." << std::endl;
#endif
          break;
        }
      }

      OutputEDep.LepExitTopology =
          (OutputEDep.LepExit ? 1 : 0) *
          (1 * OutputEDep.LepExitBack + 2 * OutputEDep.LepExitFront +
           3 * OutputEDep.LepExitYLow + 4 * OutputEDep.LepExitYHigh +
           5 * OutputEDep.LepExitXLow + 6 * OutputEDep.LepExitXHigh);
    }
    ////////End exiting lepton section
    OutputEDep.LepDep_FV = Jaccumulate(rdr->LepDep, stopBox.X_fv[0],
                                       stopBox.X_fv[1], 0, kIncludeFVYZ) +
                           Jaccumulate(rdr->LepDaughterDep, stopBox.X_fv[0],
                                       stopBox.X_fv[1], 0, kIncludeFVYZ);
    OutputEDep.LepDep_veto =
        Jaccumulate(rdr->LepDep, stopBox.X_veto_left[0], stopBox.X_veto_left[1],
                    0, kIncludeWholeDet) +
        Jaccumulate(rdr->LepDaughterDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->LepDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->LepDaughterDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->LepDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeOOFVYZ) +
        Jaccumulate(rdr->LepDaughterDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeOOFVYZ);

    OutputEDep.ProtonDep_FV =
        Jaccumulate(rdr->ProtonDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeFVYZ) +
        Jaccumulate(rdr->ProtonDaughterDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeFVYZ);
    OutputEDep.ProtonDep_veto =
        Jaccumulate(rdr->ProtonDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->ProtonDaughterDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->ProtonDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->ProtonDaughterDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->ProtonDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeOOFVYZ) +
        Jaccumulate(rdr->ProtonDaughterDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeOOFVYZ);

    OutputEDep.NeutronDep_FV =
        Jaccumulate(rdr->NeutronDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeFVYZ) +
        Jaccumulate(rdr->NeutronDaughterDep, stopBox.X_fv[0], stopBox.X_fv[1],
                    0, kIncludeFVYZ);
    OutputEDep.NeutronDep_veto =
        Jaccumulate(rdr->NeutronDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->NeutronDaughterDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->NeutronDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->NeutronDaughterDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->NeutronDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeOOFVYZ) +
        Jaccumulate(rdr->NeutronDaughterDep, stopBox.X_fv[0], stopBox.X_fv[1],
                    0, kIncludeOOFVYZ);

    OutputEDep.NeutronDep_ChrgWAvgTime_FV =
        Jaccumulate(rdr->NeutronDep_ChrgWSumTime, stopBox.X_fv[0],
                    stopBox.X_fv[1], 0, kIncludeFVYZ) +
        Jaccumulate(rdr->NeutronDaughterDep_ChrgWSumTime, stopBox.X_fv[0],
                    stopBox.X_fv[1], 0, kIncludeFVYZ);
    OutputEDep.NeutronDep_ChrgWAvgTime_veto =
        Jaccumulate(rdr->NeutronDep_ChrgWSumTime, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->NeutronDaughterDep_ChrgWSumTime,
                    stopBox.X_veto_left[0], stopBox.X_veto_left[1], 0,
                    kIncludeWholeDet) +
        Jaccumulate(rdr->NeutronDep_ChrgWSumTime, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->NeutronDaughterDep_ChrgWSumTime,
                    stopBox.X_veto_right[0], stopBox.X_veto_right[1], 0,
                    kIncludeWholeDet) +
        Jaccumulate(rdr->NeutronDep_ChrgWSumTime, stopBox.X_fv[0],
                    stopBox.X_fv[1], 0, kIncludeOOFVYZ) +
        Jaccumulate(rdr->NeutronDaughterDep_ChrgWSumTime, stopBox.X_fv[0],
                    stopBox.X_fv[1], 0, kIncludeOOFVYZ);

    OutputEDep.NeutronDep_ChrgWAvgTime_FV /= OutputEDep.NeutronDep_FV;
    OutputEDep.NeutronDep_ChrgWAvgTime_veto /= OutputEDep.NeutronDep_veto;

    OutputEDep.PiCDep_FV = Jaccumulate(rdr->PiCDep, stopBox.X_fv[0],
                                       stopBox.X_fv[1], 0, kIncludeFVYZ) +
                           Jaccumulate(rdr->PiCDaughterDep, stopBox.X_fv[0],
                                       stopBox.X_fv[1], 0, kIncludeFVYZ);
    OutputEDep.PiCDep_veto =
        Jaccumulate(rdr->PiCDep, stopBox.X_veto_left[0], stopBox.X_veto_left[1],
                    0, kIncludeWholeDet) +
        Jaccumulate(rdr->PiCDaughterDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->PiCDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->PiCDaughterDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->PiCDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeOOFVYZ) +
        Jaccumulate(rdr->PiCDaughterDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeOOFVYZ);

    OutputEDep.Pi0Dep_FV = Jaccumulate(rdr->Pi0Dep, stopBox.X_fv[0],
                                       stopBox.X_fv[1], 0, kIncludeFVYZ) +
                           Jaccumulate(rdr->Pi0DaughterDep, stopBox.X_fv[0],
                                       stopBox.X_fv[1], 0, kIncludeFVYZ);
    OutputEDep.Pi0Dep_veto =
        Jaccumulate(rdr->Pi0Dep, stopBox.X_veto_left[0], stopBox.X_veto_left[1],
                    0, kIncludeWholeDet) +
        Jaccumulate(rdr->Pi0DaughterDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->Pi0Dep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->Pi0DaughterDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->Pi0Dep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeOOFVYZ) +
        Jaccumulate(rdr->Pi0DaughterDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeOOFVYZ);

    OutputEDep.OtherDep_FV = Jaccumulate(rdr->OtherDep, stopBox.X_fv[0],
                                         stopBox.X_fv[1], 0, kIncludeFVYZ) +
                             Jaccumulate(rdr->OtherDaughterDep, stopBox.X_fv[0],
                                         stopBox.X_fv[1], 0, kIncludeFVYZ);
    OutputEDep.OtherDep_veto =
        Jaccumulate(rdr->OtherDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->OtherDaughterDep, stopBox.X_veto_left[0],
                    stopBox.X_veto_left[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->OtherDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->OtherDaughterDep, stopBox.X_veto_right[0],
                    stopBox.X_veto_right[1], 0, kIncludeWholeDet) +
        Jaccumulate(rdr->OtherDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeOOFVYZ) +
        Jaccumulate(rdr->OtherDaughterDep, stopBox.X_fv[0], stopBox.X_fv[1], 0,
                    kIncludeOOFVYZ);

    OutputEDep.TotalNonlep_Dep_FV =
        OutputEDep.ProtonDep_FV + OutputEDep.NeutronDep_FV +
        OutputEDep.PiCDep_FV + OutputEDep.Pi0Dep_FV + OutputEDep.OtherDep_FV;

    OutputEDep.TotalNonlep_Dep_veto =
        OutputEDep.ProtonDep_veto + OutputEDep.NeutronDep_veto +
        OutputEDep.PiCDep_veto + OutputEDep.Pi0Dep_veto +
        OutputEDep.OtherDep_veto;

    if (((OutputEDep.TotalNonlep_Dep_veto != 0) &&
         (!std::isnormal(OutputEDep.TotalNonlep_Dep_veto))) ||
        ((OutputEDep.TotalNonlep_Dep_FV != 0) &&
         (!std::isnormal(OutputEDep.TotalNonlep_Dep_FV)))) {
      std::cout << "\n[INFO] XBins: " << stopBox.X_fv[0] << " -- "
                << stopBox.X_fv[1] << ", Veto left: " << stopBox.X_veto_left[0]
                << " -- " << stopBox.X_veto_left[1]
                << ", Veto right: " << stopBox.X_veto_right[0] << " -- "
                << stopBox.X_veto_right[1] << std::endl;

      std::cout << "FV -- Total: " << OutputEDep.TotalNonlep_Dep_FV
                << ", Proton: " << OutputEDep.ProtonDep_FV
                << ", Neutron: " << OutputEDep.NeutronDep_FV
                << ", PiC: " << OutputEDep.PiCDep_FV
                << ", Pi0: " << OutputEDep.Pi0Dep_FV
                << ", Other: " << OutputEDep.OtherDep_FV << std::endl;

      std::cout << "Veto -- Total: " << OutputEDep.TotalNonlep_Dep_veto
                << ", Proton: " << OutputEDep.ProtonDep_veto
                << ", Neutron: " << OutputEDep.NeutronDep_veto
                << ", PiC: " << OutputEDep.PiCDep_veto
                << ", Pi0: " << OutputEDep.Pi0Dep_veto
                << ", Other: " << OutputEDep.OtherDep_veto << std::endl
                << std::endl;
    }

    OutputEDep.IsNumu = (abs(rdr->nu_PDG) == 14);
    OutputEDep.IsAntinu = (rdr->nu_PDG < 0);
    OutputEDep.IsCC = !(abs(rdr->nu_PDG) - abs(rdr->PrimaryLepPDG) - 1);
    OutputEDep.Is0Pi =
        ((rdr->NPiC + rdr->NPi0 + rdr->NGamma + rdr->NOther) == 0);
    OutputEDep.Is1PiC =
        ((rdr->NGamma + rdr->NOther + rdr->NPi0) == 0) && (rdr->NPiC == 1);
    OutputEDep.Is1Pi0 =
        ((rdr->NGamma + rdr->NOther + rdr->NPiC) == 0) && (rdr->NPi0 == 1);
    OutputEDep.Is1Pi = OutputEDep.Is1PiC || OutputEDep.Is1Pi0;
    OutputEDep.IsNPi =
        ((rdr->NGamma + rdr->NOther) == 0) && ((rdr->NPiC + rdr->NPi0) > 1);
    OutputEDep.IsOther = (rdr->NGamma + rdr->NOther);
    OutputEDep.Topology =
        (OutputEDep.IsCC ? 1 : -1) *
        (1 * OutputEDep.Is0Pi + 2 * OutputEDep.Is1PiC + 3 * OutputEDep.Is1Pi0 +
         4 * OutputEDep.IsNPi + 5 * OutputEDep.IsOther);

    OutputEDep.HadrShowerContainedInFV =
        (OutputEDep.TotalNonlep_Dep_veto < VetoThreshold);
    OutputEDep.PrimaryLeptonContainedInFV =
        (OutputEDep.LepDep_veto < VetoThreshold);

    OutputTree->Fill();
  }

  std::cout << "\r[INFO]: Read " << rdr->GetEntries() << " entries."
            << std::endl;

  outfile->Write();
  outfile->Close();
}
