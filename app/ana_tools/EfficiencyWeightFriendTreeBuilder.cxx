#include "DepositsSummaryTreeReader.hxx"
#include "StopConfigTreeReader.hxx"
#include "SelectionSummaryTreeReader.hxx"

#include "StringParserUtility.hxx"
#include "ROOTUtility.hxx"

#include "TMath.h"
#include "TTree.h"

#include <string>
#include <vector>

std::string InputDepostisSummaryFile;
std::string OutputFile;
std::string EffCorrectorFile;

bool SelMuExit = false;
double SelMuExitKE = 0;  // GeV

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "\t-i <stopprocessor.root>     : TChain descriptor for"
         " input tree. \n"
         "\t                              Should be the output of RunSelection."
         "\t-o <outputfile.root>        : Output file to write "
         "efficiency friend tree to.\n"
         "\t-E <effcorrector.root>      : File containing efficiency "
         "histograms.\n"
         "\t-m <muon exit KE>           : Muon exit threshold KE in MeV.\n"
         "\t-M <eff correction mode>    : Kinematics to use to determine the "
         "hadronic shower selection\n"
         "\t                              efficiency correction. {Default = 3}"
         "\n"
         "\t                              1: True hadronic energy, absolute "
         "off axis position.\n"
         "\t                              2: True hadronic energy, position "
         "within the detector.\n"
         "\t                              3: Visible hadronic energy, position "
         "within the detector.\n"
      << std::endl;
}

enum EffCorrModeEnum { kEHadrAbsPos = 1, kEHadrDetPos, kEHadrVisDetPos };

EffCorrModeEnum EffCorrMode = kEHadrVisDetPos;

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      InputDepostisSummaryFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      OutputFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-E") {
      EffCorrectorFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-m") {
      SelMuExitKE = str2T<double>(argv[++opt]) * 1E-3;
      SelMuExit = true;
    } else if (std::string(argv[opt]) == "-M") {
      EffCorrMode = static_cast<EffCorrModeEnum>(str2T<int>(argv[++opt]));
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  StopConfig csRdr("StopConfigTree", InputDepostisSummaryFile);
  SelectionSummary ssRdr("SelectionSummaryTree", InputDepostisSummaryFile);

  std::vector<BoundingBox> StopActiveRegions = csRdr.GetStopBoundingBoxes(false);

  TH2 *MuonKinematics_seleff =
    GetHistogram<TH2>(EffCorrectorFile, "MuonKinematics_eff");
  TH2 *ShowerKinematics_seleff = nullptr;

  switch (EffCorrMode) {
    case kEHadrAbsPos:{
      ShowerKinematics_seleff =
        GetHistogram<TH2>(EffCorrectorFile, "Ehadr_FV_abspos_seleff");
      break;
    }
    case kEHadrDetPos:{
      ShowerKinematics_seleff =
        GetHistogram<TH2>(EffCorrectorFile, "Ehadr_FV_detpos_seleff");
      break;
    }
    case kEHadrVisDetPos:{
      ShowerKinematics_seleff =
        GetHistogram<TH2>(EffCorrectorFile, "EVisHadr_FV_detpos_seleff");
      break;
    }
  }

  DepositsSummary edr("DepositsSummaryTree", InputDepostisSummaryFile);

  TFile *of = CheckOpenFile(OutputFile, "RECREATE");

  size_t loud_every = edr.GetEntries() / 10;
  Long64_t NEntries = edr.GetEntries();

  TTree *friendtree = new TTree("EffWeights", "");

  double effweight, mueffweight, hadreffweight;
  friendtree->Branch("MuEffWeight", &mueffweight, "MuEffWeight/D");
  friendtree->Branch("HadrEffWeight", &hadreffweight, "HadrEffWeight/D");
  friendtree->Branch("EffWeight", &effweight, "EffWeight/D");

  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    edr.GetEntry(e_it);

    mueffweight = 0;
    hadreffweight = 0;
    effweight = 0;

    if ((edr.stop < 0) || (edr.PrimaryLepPDG != 13)) {
      std::cout << "[ERROR]: Currently EfficiencyWeightFriendTreeBuilder "
                   "expected to receive pre-selected events."
                << std::endl;
      exit(1);
    }

    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << edr.vtx[0] << ", " << edr.vtx[1] << ", " << edr.vtx[2]
                << "}, Enu: " << edr.nu_4mom[3]
                << " ). Stop weight = " << edr.stop_weight << std::endl;
    }

    if (SelMuExit) {

      TVector3 TrStr(edr.vtx[0], edr.vtx[1], edr.vtx[2]);
      TVector3 TrDir(edr.PrimaryLep_4mom[0], edr.PrimaryLep_4mom[1],
                     edr.PrimaryLep_4mom[2]);
      TrDir = TrDir.Unit();

      double MuToWall = CalculateToWall(StopActiveRegions[edr.stop],
        TrStr, TrDir) * 1E-2;

      Int_t mu_eff_x_bin =
          MuonKinematics_seleff->GetXaxis()->FindFixBin(edr.PrimaryLep_4mom[3]);
      Int_t mu_eff_y_bin =
          MuonKinematics_seleff->GetYaxis()->FindFixBin(MuToWall);

      if (MuonKinematics_seleff->GetBinContent(mu_eff_x_bin, mu_eff_y_bin)) {
        mueffweight =
            1.0 /
            MuonKinematics_seleff->GetBinContent(mu_eff_x_bin, mu_eff_y_bin);
      } else {
        mueffweight = 1;
      }
    }


    Int_t hadr_eff_x_bin = 0;
    Int_t hadr_eff_y_bin = 0;
    switch (EffCorrMode) {
      case kEHadrAbsPos:{
        hadr_eff_x_bin = ShowerKinematics_seleff->GetXaxis()->FindFixBin(
            edr.ERecProxy_True - edr.PrimaryLep_4mom[3]);
        hadr_eff_y_bin =
            ShowerKinematics_seleff->GetYaxis()->FindFixBin(edr.vtx[0]);
        break;
      }
    case kEHadrDetPos:{
      hadr_eff_x_bin = ShowerKinematics_seleff->GetXaxis()->FindFixBin(
          edr.ERecProxy_True - edr.PrimaryLep_4mom[3]);
      hadr_eff_y_bin =
          ShowerKinematics_seleff->GetYaxis()->FindFixBin(edr.vtxInDetX);
      break;
      }
    case kEHadrVisDetPos:{
      hadr_eff_x_bin = ShowerKinematics_seleff->GetXaxis()->FindFixBin(
          edr.TotalNonlep_Dep_FV + edr.TotalNonlep_Dep_veto);
      hadr_eff_y_bin =
          ShowerKinematics_seleff->GetYaxis()->FindFixBin(edr.vtxInDetX);
      break;
      }
    }

    if (ShowerKinematics_seleff->GetBinContent(hadr_eff_x_bin,
                                              hadr_eff_y_bin)) {
      hadreffweight = 1.0 /
                      ShowerKinematics_seleff->GetBinContent(hadr_eff_x_bin,
                                                            hadr_eff_y_bin);
    } else {
      hadreffweight = 1;
    }

    effweight = (SelMuExit ? mueffweight : 1) * hadreffweight;

    friendtree->Fill();
  }

  of->Write();
  of->Close();
}
