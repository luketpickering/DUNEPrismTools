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
std::vector<double> fvgap;

TVector3 BeamPos{0, 58.1, -575.0};

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-i <fulldetprocess.root>             : TChain descriptor for"
               " input tree. \n"
               "\t-o <outputfile.root>                 : Output file to write "
               "friend tree to.\n"
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

  TChain *config_in = new TChain("stopConfigTree", "");
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

  TFile *of = new TFile(oupfile.c_str(), "RECREATE");

  TH2D *MuonKinematics_all =
      new TH2D("MuonKinematics_all", ";#it{E}_{#mu};ToWall_{#mu};Count", 100, 0,
               20, 60, 0, 6);
  TH2D *MuonKinematics_musel =
      new TH2D("MuonKinematics_musel", ";#it{E}_{#mu};ToWall_{#mu};Count", 100,
               0, 20, 60, 0, 6);
  TH2D *MuonKinematics_muhadrsel =
      new TH2D("MuonKinematics_selected", ";#it{E}_{#mu};ToWall_{#mu};Count",
               100, 0, 20, 60, 0, 6);

  TH2D *MuonKinematics_seleff =
      new TH2D("MuonKinematics_eff", ";#it{E}_{#mu};ToWall_{#mu};#epsilon", 100,
               0, 20, 60, 0, 6);

  TH2D *EDephadr_all =
      new TH2D("EDephadr_all", ";#it{EDep}_{Hadr};ToWall^{#prime}_{#mu};Count",
               100, 0, 20, 60, 0, 6);
  TH2D *EDephadr_hadrsel = new TH2D(
      "EDephadr_hadrsel", ";#it{EDep}_{Hadr};ToWall^{#prime}_{#mu};Count", 100,
      0, 20, 60, 0, 6);
  TH2D *EDephadr_muhadrsel = new TH2D(
      "EDephadr_selected", ";#it{EDep}_{Hadr};ToWall^{#prime}_{#mu};Count", 100,
      0, 20, 60, 0, 6);

  TH2D *EDephadr_seleff = new TH2D(
      "EDephadr_eff", ";#it{EDep}_{Hadr};ToWall^{#prime}_{#mu};#epsilon", 100,
      0, 20, 60, 0, 6);

  TH2D *EDephadr_true_all =
      new TH2D("EDephadr_true_all",
               ";#it{EDep}_{Hadr};ToWall_{Charged hadrons true};Count", 100, 0,
               20, 60, 0, 6);
  TH2D *EDephadr_true_hadrsel =
      new TH2D("EDephadr_true_hadrsel",
               ";#it{EDep}_{Hadr};ToWall_{Charged hadrons true};Count", 100, 0,
               20, 60, 0, 6);
  TH2D *EDephadr_true_muhadrsel =
      new TH2D("EDephadr_true_selected",
               ";#it{EDep}_{Hadr};ToWall_{Charged hadrons true};Count", 100, 0,
               20, 60, 0, 6);

  TH2D *EDephadr_true_seleff =
      new TH2D("EDephadr_true_eff",
               ";#it{EDep}_{Hadr};ToWall_{Charged hadrons true};#epsilon", 100,
               0, 20, 60, 0, 6);

  TH2D *EDephadr_offset_all = new TH2D(
      "EDephadr_offset_all", ";#it{EDep}_{Hadr};Off axis position (cm);Count",
      100, 0, 20, 80, -3750, 150);
  TH2D *EDephadr_offset_muhadrsel =
      new TH2D("EDephadr_offset_muhadrsel",
               ";#it{EDep}_{Hadr};Off axis position (cm);Count", 100, 0, 20, 80,
               -3750, 150);

  TH2D *EDephadr_offset_seleff =
      new TH2D("EDephadr_offset_seleff",
               ";#it{EDep}_{Hadr};Off axis position (cm);#epsilon", 100, 0, 20,
               80, -3750, 150);

  TH2D *EDephadr_FV_offset_all =
      new TH2D("EDephadr_FV_offset_all",
               ";#it{EDep}_{Hadr};Off axis position (cm);Count", 100, 0, 20, 80,
               -3750, 150);
  TH2D *EDephadr_FV_offset_muhadrsel =
      new TH2D("EDephadr_FV_offset_muhadrsel",
               ";#it{EDep}_{Hadr};Off axis position (cm);Count", 100, 0, 20, 80,
               -3750, 150);

  TH2D *EDephadr_FV_offset_seleff =
      new TH2D("EDephadr_FV_offset_seleff",
               ";#it{EDep}_{Hadr};Off axis position (cm);#epsilon", 100, 0, 20,
               80, -3750, 150);

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

    TVector3 TrDir_rot = TrDir;
    TrDir_rot.Rotate(TMath::Pi(), (TrStr - BeamPos).Unit());

    double ToWall_hadr =
        CalculateToWall(BBs[edr.stop], TrStr, TrDir_rot) * 1E-2;

    EDephadr_all->Fill(edr.TotalNonlep_Dep_FV + edr.TotalNonlep_Dep_veto,
                       ToWall_hadr);
    EDephadr_hadrsel->Fill(edr.TotalNonlep_Dep_FV + edr.TotalNonlep_Dep_veto,
                           ToWall_hadr, edr.HadrShowerContainedInFV);
    EDephadr_muhadrsel->Fill(
        edr.TotalNonlep_Dep_FV + edr.TotalNonlep_Dep_veto, ToWall_hadr,
        edr.LepExit_AboveThresh * edr.HadrShowerContainedInFV);

    TVector3 HadrShower_dir(0, 0, 0);
    for (size_t pt_it = 0; pt_it < edr.GetNPassthroughParts(); ++pt_it) {
      std::pair<Int_t, Double_t *> part = edr.GetPassthroughPart(pt_it);

      if ((abs(part.first) > 9999) || (abs(part.first) == 14) ||
          (abs(part.first) == 12) || (abs(part.first) == 2112)) {
        continue;
      }

      HadrShower_dir +=
          TVector3(part.second[0], part.second[1], part.second[2]);
    }

    double ToWall_hadr_true =
        CalculateToWall(BBs[edr.stop], TrStr, HadrShower_dir) * 1E-2;

    EDephadr_true_all->Fill(edr.TotalNonlep_Dep_FV + edr.TotalNonlep_Dep_veto,
                            ToWall_hadr_true);
    EDephadr_true_hadrsel->Fill(
        edr.TotalNonlep_Dep_FV + edr.TotalNonlep_Dep_veto, ToWall_hadr_true,
        edr.HadrShowerContainedInFV);
    EDephadr_true_muhadrsel->Fill(
        edr.TotalNonlep_Dep_FV + edr.TotalNonlep_Dep_veto, ToWall_hadr_true,
        edr.LepExit_AboveThresh * edr.HadrShowerContainedInFV);

    EDephadr_offset_all->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                                  edr.TotalNonlep_Dep_veto,
                              edr.vtx[0]);
    EDephadr_offset_muhadrsel->Fill(
        edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
            edr.TotalNonlep_Dep_veto,
        edr.vtx[0], edr.LepExit_AboveThresh * edr.HadrShowerContainedInFV);

    EDephadr_FV_offset_all->Fill(
        edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV, edr.vtx[0]);
    EDephadr_FV_offset_muhadrsel->Fill(
        edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV, edr.vtx[0],
        edr.LepExit_AboveThresh * edr.HadrShowerContainedInFV);
  }

  MuonKinematics_seleff->Divide(MuonKinematics_musel, MuonKinematics_all);
  EDephadr_seleff->Divide(EDephadr_hadrsel, EDephadr_all);
  EDephadr_true_seleff->Divide(EDephadr_true_hadrsel, EDephadr_true_all);
  EDephadr_offset_seleff->Divide(EDephadr_offset_muhadrsel,
                                 EDephadr_offset_all);
  EDephadr_FV_offset_seleff->Divide(EDephadr_FV_offset_muhadrsel,
                                 EDephadr_FV_offset_all);

  TTree *friendtree = new TTree("EffWeights", "");

  double effweight, mueffweight, hadreffweight, MuToWall, HadrToWall,
      hadreffweight_true, HadrToWall_true, dumbeffweight, dumbeffweight_FV;
  double HadrShower_dir_arr[3], RotMu_dir_arr[3], Mu_dir_arr[3];
  friendtree->Branch("MuEffWeight", &mueffweight, "MuEffWeight/D");
  friendtree->Branch("HadrEffWeight", &hadreffweight, "HadrEffWeight/D");
  friendtree->Branch("HadrEffWeight_true", &hadreffweight_true,
                     "HadrEffWeight_true/D");
  friendtree->Branch("EffWeight", &effweight, "EffWeight/D");
  friendtree->Branch("MuToWall", &MuToWall, "MuToWall/D");
  friendtree->Branch("HadrToWall", &HadrToWall, "HadrToWall/D");
  friendtree->Branch("HadrToWall_true", &HadrToWall_true, "HadrToWall_true/D");
  friendtree->Branch("HadrShower_dir", &HadrShower_dir_arr,
                     "HadrShower_dir[3]/D");
  friendtree->Branch("RotMu_dir", &RotMu_dir_arr, "RotMu_dir[3]/D");
  friendtree->Branch("Mu_dir", &Mu_dir_arr, "Mu_dir[3]/D");
  friendtree->Branch("DumbEffWeight", &dumbeffweight, "DumbEffWeight/D");
  friendtree->Branch("DumbEffWeight_FV", &dumbeffweight_FV,
                     "DumbEffWeight_FV/D");

  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    edr.GetEntry(e_it);
    effweight = 0;
    mueffweight = 0;
    hadreffweight = 0;
    MuToWall = 0;
    HadrToWall = 0;
    std::fill_n(Mu_dir_arr, 3, 0);
    std::fill_n(RotMu_dir_arr, 3, 0);
    std::fill_n(HadrShower_dir_arr, 3, 0);
    HadrToWall_true = 0;
    hadreffweight_true = 0;

    if ((edr.stop == -1) || (edr.PrimaryLepPDG != 13)) {
      friendtree->Fill();
      continue;
    }

    TVector3 TrStr(edr.vtx[0], edr.vtx[1], edr.vtx[2]);
    TVector3 TrDir(edr.PrimaryLep_4mom[0], edr.PrimaryLep_4mom[1],
                   edr.PrimaryLep_4mom[2]);
    TrDir = TrDir.Unit();
    Mu_dir_arr[0] = TrDir[0];
    Mu_dir_arr[1] = TrDir[1];
    Mu_dir_arr[2] = TrDir[2];

    MuToWall = CalculateToWall(BBs[edr.stop], TrStr, TrDir) * 1E-2;

    TVector3 TrDir_rot = TrDir;
    TrDir_rot.Rotate(TMath::Pi(), (TrStr - BeamPos).Unit());

    RotMu_dir_arr[0] = TrDir_rot[0];
    RotMu_dir_arr[1] = TrDir_rot[1];
    RotMu_dir_arr[2] = TrDir_rot[2];

    HadrToWall = CalculateToWall(BBs[edr.stop], TrStr, TrDir_rot) * 1E-2;

    TVector3 HadrShower_dir(0, 0, 0);
    for (size_t pt_it = 0; pt_it < edr.GetNPassthroughParts(); ++pt_it) {
      std::pair<Int_t, Double_t *> part = edr.GetPassthroughPart(pt_it);

      if ((abs(part.first) > 9999) || (abs(part.first) == 14) ||
          (abs(part.first) == 12) || (abs(part.first) == 2112)) {
        continue;
      }

      HadrShower_dir +=
          TVector3(part.second[0], part.second[1], part.second[2]);
    }
    HadrShower_dir = HadrShower_dir.Unit();
    HadrShower_dir_arr[0] = HadrShower_dir[0];
    HadrShower_dir_arr[1] = HadrShower_dir[1];
    HadrShower_dir_arr[2] = HadrShower_dir[2];

    HadrToWall_true =
        CalculateToWall(BBs[edr.stop], TrStr, HadrShower_dir) * 1E-2;

    Int_t twmu_bin = MuonKinematics_seleff->GetYaxis()->FindFixBin(MuToWall);
    Int_t Emu_bin =
        MuonKinematics_seleff->GetXaxis()->FindFixBin(edr.PrimaryLep_4mom[3]);

    Int_t twhadr_bin = EDephadr_seleff->GetYaxis()->FindFixBin(HadrToWall);
    Int_t Ehadr_bin = EDephadr_seleff->GetXaxis()->FindFixBin(
        edr.TotalNonlep_Dep_FV + edr.TotalNonlep_Dep_veto);

    Int_t twhadr_true_bin =
        EDephadr_true_seleff->GetYaxis()->FindFixBin(HadrToWall_true);
    Int_t Ehadr_true_bin = EDephadr_true_seleff->GetXaxis()->FindFixBin(
        edr.TotalNonlep_Dep_FV + edr.TotalNonlep_Dep_veto);

    Int_t dumbeff_x_bin =
        EDephadr_offset_seleff->GetYaxis()->FindFixBin(edr.vtx[0]);
    Int_t dumbeff_edep_bin = EDephadr_offset_seleff->GetXaxis()->FindFixBin(
        edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
        edr.TotalNonlep_Dep_veto);

    Int_t dumbeff_fv_x_bin =
        EDephadr_FV_offset_seleff->GetYaxis()->FindFixBin(edr.vtx[0]);
    Int_t dumbeff_fv_edep_bin =
        EDephadr_FV_offset_seleff->GetXaxis()->FindFixBin(
            edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV);

    if (MuonKinematics_seleff->GetBinContent(Emu_bin, twmu_bin)) {
      mueffweight =
          1.0 / MuonKinematics_seleff->GetBinContent(Emu_bin, twmu_bin);
    } else {
      mueffweight = 1;
    }

    if (EDephadr_seleff->GetBinContent(Ehadr_bin, twhadr_bin)) {
      hadreffweight =
          1.0 / EDephadr_seleff->GetBinContent(Ehadr_bin, twhadr_bin);
    } else {
      hadreffweight = 1;
    }

    if (EDephadr_true_seleff->GetBinContent(Ehadr_true_bin, twhadr_true_bin)) {
      hadreffweight_true =
          1.0 /
          EDephadr_true_seleff->GetBinContent(Ehadr_true_bin, twhadr_true_bin);
    } else {
      hadreffweight_true = 1;
    }

    if (EDephadr_offset_seleff->GetBinContent(dumbeff_edep_bin,
                                              dumbeff_x_bin)) {
      dumbeffweight = 1.0 /
                      EDephadr_offset_seleff->GetBinContent(dumbeff_edep_bin,
                                                            dumbeff_x_bin);
    } else {
      dumbeffweight = 1;
    }

    if (EDephadr_FV_offset_seleff->GetBinContent(dumbeff_fv_edep_bin,
                                                 dumbeff_fv_x_bin)) {
      dumbeffweight_FV = 1.0 /
                         EDephadr_FV_offset_seleff->GetBinContent(
                             dumbeff_fv_edep_bin, dumbeff_fv_x_bin);
    } else {
      dumbeffweight_FV = 1;
    }

    effweight = mueffweight * hadreffweight;

    friendtree->Fill();
  }

  of->Write();
  of->Close();
}
