#include "BoundingBox.hxx"
#include "DepositsSummaryTreeReader.hxx"
#include "StopConfigTreeReader.hxx"

#include "StringParserUtility.hxx"
#include "ROOTUtility.hxx"
#include "GetUsage.hxx"

#include "TMath.h"
#include "TTree.h"
#include "TH2D.h"
#include "TFile.h"

#include <string>
#include <vector>

std::string InputDepostisSummaryFile;
std::string OutputFile;

double HadrVeto = 10E-3;     // GeV
double SelMuExitKE = 50E-3;  // GeV

std::vector<double> ERecBinning;
bool BuildMissingProtonEFakeData = false;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0] << "\n"
            << GetUsageText(argv[0], "ana_tools") << std::endl;
}

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
    } else if (std::string(argv[opt]) == "-v") {
      HadrVeto = str2T<double>(argv[++opt]) * 1E-3;
    } else if (std::string(argv[opt]) == "-m") {
      SelMuExitKE = str2T<double>(argv[++opt]) * 1E-3;
    } else if (std::string(argv[opt]) == "-b") {
      ERecBinning = BuildBinEdges(argv[++opt]);
    } else if (std::string(argv[opt]) == "-FDproton") {
      BuildMissingProtonEFakeData = true;
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

  if (!ERecBinning.size()) {
    int argc_dum = 3;
    char const *argv_dum[] = {"", "-b", "0_10:0.1"};
    handleOpts(argc_dum, argv_dum);
  }

  StopConfig csRdr(InputDepostisSummaryFile);

  std::vector<BoundingBox> StopActiveRegions = csRdr.GetStopBoundingBoxes(false);

  Double_t MaxStopWidth = -std::numeric_limits<double>::max();
  Double_t MaxToWall = -std::numeric_limits<double>::max();
  Double_t MaxAbsX = -std::numeric_limits<double>::max();
  Double_t MinAbsX = std::numeric_limits<double>::max();
  for(BoundingBox const & bb : StopActiveRegions){
    MaxStopWidth = std::max(MaxStopWidth, bb.Max[0] - bb.Min[0]);
    MaxToWall = std::max(MaxToWall,
      ( (TVector3(bb.Max.data()) -
         TVector3(bb.Min.data()) ).Mag()));
     MaxAbsX = std::max(MaxAbsX,bb.Max[0]);
     MinAbsX = std::min(MinAbsX,bb.Min[0]);
  }

  std::cout << "[INFO]: Max DetXRange = " << MaxStopWidth
            << ", Max ToWall = " << (MaxToWall * 1E-2) << std::endl;

  DepositsSummary edr(InputDepostisSummaryFile);

  Double_t ProtonFakeDataWeight = 1;
  if (BuildMissingProtonEFakeData) {
    // Should probably check that this exists, but leaving hardcoded and
    // presumptive for now as CheckTTreeHasBranch doesn't seem to work on
    // pre-friended trees.
    edr.tree->SetBranchAddress("EnuTp", &ProtonFakeDataWeight);
  }

  TFile *of = CheckOpenFile(OutputFile, "RECREATE");

  TH2D *MuonKinematics_all =
      new TH2D("MuonKinematics_all", ";#it{E}_{#mu};ToWall_{#mu} (m);Count",
               (ERecBinning.size()-1),ERecBinning.data(), 60, 0,
               MaxToWall * 1E-2);
  TH2D *MuonKinematics_musel =
      new TH2D("MuonKinematics_musel", ";#it{E}_{#mu};ToWall_{#mu} (m);Count",
               (ERecBinning.size()-1),ERecBinning.data(), 60, 0,
               MaxToWall * 1E-2);
  TH2D *MuonKinematics_seleff =
      new TH2D("MuonKinematics_eff", ";#it{E}_{#mu};ToWall_{#mu} (m);#epsilon",
               (ERecBinning.size()-1),ERecBinning.data(), 60, 0,
               MaxToWall * 1E-2);

  TH2D *Ehadr_FV_detpos_all = new TH2D(
      "Ehadr_FV_detpos_all", ";#it{E}_{Hadr};Detector X position (cm);Count",
      (ERecBinning.size()-1),ERecBinning.data(), 100, -(MaxStopWidth / 2.0),
      (MaxStopWidth / 2.0));
  TH2D *Ehadr_FV_detpos_hadrsel =
      new TH2D("Ehadr_FV_detpos_hadrsel",
               ";#it{E}_{Hadr};Detector X position (cm);Count",
               (ERecBinning.size()-1),ERecBinning.data(), 100,
               -(MaxStopWidth / 2.0), (MaxStopWidth / 2.0));
  TH2D *Ehadr_FV_detpos_seleff =
      new TH2D("Ehadr_FV_detpos_seleff",
               ";#it{E}_{Hadr};Detector X position (cm);#epsilon",
               (ERecBinning.size()-1),ERecBinning.data(),
               100, -(MaxStopWidth / 2.0), (MaxStopWidth / 2.0));

   TH2D *EVisHadr_FV_detpos_all = new TH2D(
       "EVisHadr_FV_detpos_all",
       ";#it{E}_{Vis,Hadr};Detector X position (cm);Count",
       (ERecBinning.size()-1),ERecBinning.data(), 100, -(MaxStopWidth / 2.0),
       (MaxStopWidth / 2.0));
   TH2D *EVisHadr_FV_detpos_hadrsel =
       new TH2D("EVisHadr_FV_detpos_hadrsel",
                ";#it{E}_{Vis,Hadr};Detector X position (cm);Count",
                (ERecBinning.size()-1),ERecBinning.data(), 100,
                -(MaxStopWidth / 2.0), (MaxStopWidth / 2.0));
   TH2D *EVisHadr_FV_detpos_seleff =
       new TH2D("EVisHadr_FV_detpos_seleff",
                ";#it{E}_{Vis,Hadr};Detector X position (cm);#epsilon",
                (ERecBinning.size()-1),ERecBinning.data(),
                100, -(MaxStopWidth / 2.0), (MaxStopWidth / 2.0));

    TH2D *ETrueNonNeutron_FV_detpos_all = new TH2D(
        "ETrueNonNeutron_FV_detpos_all",
        ";#it{E}_{Vis,Hadr};Detector X position (cm);Count",
        (ERecBinning.size()-1),ERecBinning.data(), 100, -(MaxStopWidth / 2.0),
        (MaxStopWidth / 2.0));
    TH2D *ETrueNonNeutron_FV_detpos_hadrsel =
        new TH2D("ETrueNonNeutron_FV_detpos_hadrsel",
                 ";#it{E}_{Vis,Hadr};Detector X position (cm);Count",
                 (ERecBinning.size()-1),ERecBinning.data(), 100,
                 -(MaxStopWidth / 2.0), (MaxStopWidth / 2.0));
    TH2D *ETrueNonNeutron_FV_detpos_seleff =
        new TH2D("ETrueNonNeutron_FV_detpos_seleff",
                 ";#it{E}_{Vis,Hadr};Detector X position (cm);#epsilon",
                 (ERecBinning.size()-1),ERecBinning.data(),
                 100, -(MaxStopWidth / 2.0), (MaxStopWidth / 2.0));


   TH2D *Ehadr_FV_abspos_all = new TH2D(
       "Ehadr_FV_abspos_all", ";#it{E}_{Hadr};Detector X position (cm);Count",
       (ERecBinning.size()-1),ERecBinning.data(), 100, MinAbsX, MaxAbsX);
   TH2D *Ehadr_FV_abspos_hadrsel =
       new TH2D("Ehadr_FV_abspos_hadrsel",
                ";#it{E}_{Hadr};Detector X position (cm);Count",
                (ERecBinning.size()-1),ERecBinning.data(), 100,
                MinAbsX, MaxAbsX);
   TH2D *Ehadr_FV_abspos_seleff =
       new TH2D("Ehadr_FV_abspos_seleff",
                ";#it{E}_{Hadr};Detector X position (cm);#epsilon",
                (ERecBinning.size()-1),ERecBinning.data(),
                100, MinAbsX, MaxAbsX);

  size_t loud_every = edr.GetEntries() / 10;
  Long64_t NEntries = edr.GetEntries();

  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    edr.GetEntry(e_it);
    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << edr.vtx[0] << ", " << edr.vtx[1] << ", " << edr.vtx[2]
                << "}, Enu: " << edr.nu_4mom[3] << " )" << std::endl;

      if (BuildMissingProtonEFakeData) {
        std::cout << "\tProton fake data weight = " << ProtonFakeDataWeight
                  << std::endl;
      }
    }

    if ((edr.stop == -1) || (edr.PrimaryLepPDG != 13)) {
      continue;
    }

    if (BuildMissingProtonEFakeData) {
      edr.TotalNonlep_Dep_veto -= edr.ProtonDep_veto * 0.2;
      edr.ERecProxy_True -= edr.EKinProton_True * 0.2;
    }

    TVector3 TrStr(edr.vtx[0], edr.vtx[1], edr.vtx[2]);
    TVector3 TrDir(edr.PrimaryLep_4mom[0], edr.PrimaryLep_4mom[1],
                   edr.PrimaryLep_4mom[2]);
    TrDir = TrDir.Unit();

    double ToWall = CalculateToWall(StopActiveRegions[edr.stop],
      TrStr, TrDir) * 1E-2;

    MuonKinematics_all->Fill(edr.GetProjection(DepositsSummary::kEFSLep_True),
      ToWall,ProtonFakeDataWeight);
    MuonKinematics_musel->Fill(
        edr.GetProjection(DepositsSummary::kEFSLep_True), ToWall,
        ProtonFakeDataWeight * double(edr.LepExitKE > SelMuExitKE));

    Ehadr_FV_detpos_all->Fill(edr.GetProjection(DepositsSummary::kEHadr_True),
                              edr.vtxInDetX, ProtonFakeDataWeight);
    Ehadr_FV_detpos_hadrsel->Fill(
        edr.GetProjection(DepositsSummary::kEHadr_True), edr.vtxInDetX,
        ProtonFakeDataWeight * double(edr.TotalNonlep_Dep_veto < HadrVeto));

    EVisHadr_FV_detpos_all->Fill(
      edr.GetProjection(DepositsSummary::kEHadr_vis), edr.vtxInDetX,
      ProtonFakeDataWeight);
    EVisHadr_FV_detpos_hadrsel->Fill(
      edr.GetProjection(DepositsSummary::kEHadr_vis), edr.vtxInDetX,
      ProtonFakeDataWeight * double(edr.TotalNonlep_Dep_veto < HadrVeto));

    Ehadr_FV_abspos_all->Fill(edr.GetProjection(DepositsSummary::kEHadr_True),
      edr.vtx[0], ProtonFakeDataWeight);
    Ehadr_FV_abspos_hadrsel->Fill(edr.GetProjection(DepositsSummary::kEHadr_True),
      edr.vtx[0],
      ProtonFakeDataWeight * double(edr.TotalNonlep_Dep_veto < HadrVeto));


    ETrueNonNeutron_FV_detpos_all->Fill(edr.GetProjection(DepositsSummary::kENonNeutronHadr_True),
      edr.vtx[0], ProtonFakeDataWeight);
    ETrueNonNeutron_FV_detpos_hadrsel->Fill(edr.GetProjection(DepositsSummary::kENonNeutronHadr_True),
      edr.vtx[0],
      ProtonFakeDataWeight * double(edr.TotalNonlep_Dep_veto < HadrVeto));
  }

  MuonKinematics_seleff->Divide(MuonKinematics_musel, MuonKinematics_all);
  Ehadr_FV_detpos_seleff->Divide(Ehadr_FV_detpos_hadrsel, Ehadr_FV_detpos_all);
  EVisHadr_FV_detpos_seleff->Divide(EVisHadr_FV_detpos_hadrsel,
    EVisHadr_FV_detpos_all);
  Ehadr_FV_abspos_seleff->Divide(Ehadr_FV_abspos_hadrsel, Ehadr_FV_abspos_all);
  ETrueNonNeutron_FV_detpos_seleff->Divide(ETrueNonNeutron_FV_detpos_hadrsel, ETrueNonNeutron_FV_detpos_all);

  of->Write();
  of->Close();
}
