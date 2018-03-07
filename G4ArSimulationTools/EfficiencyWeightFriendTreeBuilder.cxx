#include "BoundingBox.h"
#include "EDepTreeReader.h"

#include "DetectorStop.hxx"
#include "Utils.hxx"

#include "TTree.h"

#include <string>
#include <vector>

std::vector<DetectorStop> DetectorStops;
std::vector<BoundingBox> BBs;

std::string inpfile;
std::string oupfile;
std::string runPlanCfg, runPlanName = "";
std::vector<double> fvgap;

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "\t-i <fulldetprocess.root>             : TChain descriptor for"
         " input tree. \n"
         "\t-o <outputfile.root>                 : Output file to write "
         "friend tree to.\n"
         "\t-r <RunPlan.XML>              : An XML file specifying a run "
         "plan  \n"
         "\t                                to build fluxes for. See     "
         "      \n"
         "\t                                documentation for XML "
         "structure.   \n"
         "\t-fv <x,y,z>                  : FV gap in cm, default = 50,50,50."
      << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-i") {
      inpfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      oupfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-fv") {
      fvgap = ParseToVect<double>(argv[++opt], ",");
    } else if (std::string(argv[opt]) == "-r") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");

      runPlanCfg = params[0];

      if (params.size() > 1) {
        runPlanName = params[1];
      }
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

  DetectorStops = ReadDetectorStopConfig(runPlanCfg, runPlanName);

  size_t NDets = DetectorStops.size();

  for (size_t d_it = 0; d_it < NDets; ++d_it) {
    TVector3 Min, Max;

    Min[0] = -DetectorStops[d_it].LateralOffset * 100.0 -
             DetectorStops[d_it].DetectorFiducialWidth * 100.0 - fvgap[0];
    Min[1] = -DetectorStops[d_it].DetectorFiducialHeight * 100.0 - fvgap[1];
    Min[2] = -DetectorStops[d_it].DetectorFiducialDepth * 100.0 - fvgap[2];

    Max[0] = -DetectorStops[d_it].LateralOffset * 100.0 +
             DetectorStops[d_it].DetectorFiducialWidth * 100.0 + fvgap[0];
    Max[1] = DetectorStops[d_it].DetectorFiducialHeight * 100.0 + fvgap[1];
    Max[2] = DetectorStops[d_it].DetectorFiducialDepth * 100.0 + fvgap[2];

    BBs.emplace_back(Max, Min);
  }

  EDep edr("EDeps", inpfile);

  TFile *of = new TFile(oupfile.c_str(), "RECREATE");

  TH2D *MuonKinematics_all =
      new TH2D("MuonKinematics_all", ";#it{E}_{#mu};ToWall_{#mu};Count", 100, 0,
               5, 100, 0, 4);
  TH2D *MuonKinematics_musel =
      new TH2D("MuonKinematics_musel", ";#it{E}_{#mu};ToWall_{#mu};Count", 100,
               0, 5, 100, 0, 4);
  TH2D *MuonKinematics_muhadrsel =
      new TH2D("MuonKinematics_muhadrsel", ";#it{E}_{#mu};ToWall_{#mu};Count",
               100, 0, 5, 100, 0, 4);

  size_t loud_every = edr.GetEntries() / 10;
  Long64_t NEntries = edr.GetEntries();

  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    edr.GetEntry(e_it);
    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << edr.vtx[0] << ", " << edr.vtx[1] << ", " << edr.vtx[2]
                << "}, Enu: " << edr.nu_4mom[3] << " )" << std::endl;
    }

    if ((edr.stop == -1) || (edr.PrimaryLepPDG != 13)) {
      continue;
    }

    TVector3 TrStr(edr.vtx[0], edr.vtx[1], edr.vtx[2]);
    TVector3 TrDir(edr.PrimaryLep_4mom[0], edr.PrimaryLep_4mom[1],
                   edr.PrimaryLep_4mom[2]);
    TrDir = TrDir.Unit();

    double ToWall = CalculateToWall(BBs[edr.stop], TrStr, TrDir) * 1E-2;

    MuonKinematics_all->Fill(edr.PrimaryLep_4mom[3], ToWall);
    MuonKinematics_musel->Fill(edr.PrimaryLep_4mom[3], ToWall,
                               edr.LepExit_AboveThresh);
    MuonKinematics_muhadrsel->Fill(
        edr.PrimaryLep_4mom[3], ToWall,
        edr.LepExit_AboveThresh * edr.HadrShowerContainedInFV);
  }

  of->Write();
  of->Close();
}
