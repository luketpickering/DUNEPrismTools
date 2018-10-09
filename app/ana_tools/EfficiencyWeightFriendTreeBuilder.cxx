#include "DepositsSummaryTreeReader.hxx"
#include "StopConfigTreeReader.hxx"
#include "SelectionSummaryTreeReader.hxx"

#include "EffCorrector.hxx"

#include "StringParserUtility.hxx"
#include "ROOTUtility.hxx"
#include "GetUsage.hxx"

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
  std::cout << "[USAGE]: " << argv[0] << "\n"
            << GetUsageText(argv[0], "ana_tools") << std::endl;
}

EffCorrector::ModeEnum EffCorrMode = EffCorrector::kENonNeutronHadrDetPos;

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
      EffCorrMode = static_cast<EffCorrector::ModeEnum>(str2T<int>(argv[++opt]));
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

  StopConfig csRdr(InputDepostisSummaryFile);

  EffCorrector EffCorr(EffCorrMode, EffCorrectorFile, csRdr);

  DepositsSummary edr(InputDepostisSummaryFile);

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
      mueffweight = EffCorr.GetMuonKinematicsEffWeight(edr);
    }

    hadreffweight = EffCorr.GetHadronKinematicsEffWeight(edr);

    effweight = (SelMuExit ? mueffweight : 1) * hadreffweight;

    friendtree->Fill();
  }

  of->Write();
  of->Close();
}
