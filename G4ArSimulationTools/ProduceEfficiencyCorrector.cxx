#include "BoundingBox.h"
#include "EDepTreeReader.h"

#include "DetectorStop.hxx"
#include "Utils.hxx"

#include "TMath.h"
#include "TTree.h"

#include <string>
#include <vector>

std::vector<DetectorStop> DetectorStops;
std::vector<BoundingBox> BBs;

std::string inpfile;
std::string oupfile;

double HadrVeto = 10E-3;     // GeV
double SelMuExitKE = 50E-3;  // GeV

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "\t-i <fulldetprocess.root>    : TChain descriptor for"
         " input tree. \n"
         "\t-o <outputfile.root>        : Output file to write "
         "selected tree to.\n"
         "\t-v <hadr veto threshold>    : Hadronic shower veto threshold in "
         "MeV {default: 10}.\n"
         "\t-m <muon exit KE>           : Muon exit threshold KE in MeV "
         "{default: 0}.\n"
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
    } else if (std::string(argv[opt]) == "-v") {
      HadrVeto = str2T<double>(argv[++opt]) * 1E-3;
    } else if (std::string(argv[opt]) == "-m") {
      SelMuExitKE = str2T<double>(argv[++opt]) * 1E-3;
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

  TChain *config_in = new TChain("configTree", "");
  config_in->Add(inpfile.c_str());
  Int_t NDets;
  Double_t FVGap[3];
  config_in->SetBranchAddress("NStops", &NDets);
  config_in->SetBranchAddress("FVGap", &FVGap);
  config_in->GetEntry(0);

  TChain *stopConfig_in = new TChain("stopConfigTree", "");
  stopConfig_in->Add(inpfile.c_str());
  Double_t StopMin[3], StopMax[3];
  Double_t Offset;
  stopConfig_in->SetBranchAddress("Min", &StopMin);
  stopConfig_in->SetBranchAddress("Max", &StopMax);
  stopConfig_in->SetBranchAddress("Offset", &Offset);

  for (Int_t d_it = 0; d_it < NDets; ++d_it) {
    stopConfig_in->GetEntry(d_it);
    TVector3 Min, Max;

    Min[0] = StopMin[0];
    Min[1] = StopMin[1];
    Min[2] = StopMin[2];

    Max[0] = StopMax[0];
    Max[1] = StopMax[1];
    Max[2] = StopMax[2];

    BBs.emplace_back(Max, Min);

    std::cout << "[INFO]: Read det stop " << d_it << " {" << Min[0] << ", "
              << Min[1] << ", " << Min[2] << "} -- {" << Max[1] << ", "
              << Max[0] << ", " << Max[2] << "}" << std::endl;
  }

  EDep edr("EDeps", inpfile);

  TFile *of = CheckOpenFile(oupfile, "RECREATE");

  TH2D *MuonKinematics_all =
      new TH2D("MuonKinematics_all", ";#it{E}_{#mu};ToWall_{#mu};Count", 100, 0,
               20, 60, 0, 6);
  TH2D *MuonKinematics_musel =
      new TH2D("MuonKinematics_musel", ";#it{E}_{#mu};ToWall_{#mu};Count", 100,
               0, 20, 60, 0, 6);
  TH2D *MuonKinematics_seleff =
      new TH2D("MuonKinematics_eff", ";#it{E}_{#mu};ToWall_{#mu};#epsilon", 100,
               0, 20, 60, 0, 6);

  TH2D *Ehadr_FV_detpos_all =
      new TH2D("Ehadr_FV_detpos_all",
               ";#it{EHadr}_{Hadr};Detector X position (cm);Count", 100, 0, 20,
               120, -200, 200);
  TH2D *Ehadr_FV_detpos_hadrsel =
      new TH2D("Ehadr_FV_detpos_hadrsel",
               ";#it{EHadr}_{Hadr};Detector X position (cm);Count", 100, 0, 20,
               120, -200, 200);
  TH2D *Ehadr_FV_detpos_seleff =
      new TH2D("Ehadr_FV_detpos_seleff",
               ";#it{EHadr}_{Hadr};Detector X position (cm);#epsilon", 100, 0,
               20, 120, -200, 200);

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
                               (edr.LepExitKE > SelMuExitKE));

    Ehadr_FV_detpos_all->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                              edr.vtxInDetX);
    Ehadr_FV_detpos_hadrsel->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                                  edr.vtxInDetX,
                                  (edr.TotalNonlep_Dep_veto < HadrVeto));
  }

  MuonKinematics_seleff->Divide(MuonKinematics_musel, MuonKinematics_all);
  Ehadr_FV_detpos_seleff->Divide(Ehadr_FV_detpos_hadrsel, Ehadr_FV_detpos_all);

  of->Write();
  of->Close();
}
