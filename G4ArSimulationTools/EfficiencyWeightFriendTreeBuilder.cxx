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
std::string EffCorrectorFile;

double HadrVeto = 0;  // GeV
bool SelMuExit = false;
double SelMuExitKE = 0;  // GeV

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "\t-i <fulldetprocess.root>    : TChain descriptor for"
         " input tree. \n"
         "\t-o <outputfile.root>        : Output file to write "
         "selected tree to.\n"
         "\t-E <effcorrector.root>      : File containing efficiency "
         "histograms.\n"
         "\t-v <hadr veto threshold>    : Hadronic shower veto threshold in "
         "MeV.\n"
         "\t-m <muon exit KE>           : Muon exit threshold KE in MeV.\n"
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
    }else if (std::string(argv[opt]) == "-E") {
      EffCorrectorFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-v") {
      HadrVeto = str2T<double>(argv[++opt]) * 1E-3;
    } else if (std::string(argv[opt]) == "-m") {
      SelMuExitKE = str2T<double>(argv[++opt]) * 1E-3;
      SelMuExit = true;
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
  Int_t NStops;
  Double_t FVGap[3];
  config_in->SetBranchAddress("NStops", &NStops);
  config_in->SetBranchAddress("FVGap", &FVGap);
  config_in->GetEntry(0);

  TChain *stopConfig_in = new TChain("stopConfigTree", "");
  stopConfig_in->Add(inpfile.c_str());
  Double_t StopMin[3], StopMax[3];
  Double_t Offset;
  stopConfig_in->SetBranchAddress("Min", &StopMin);
  stopConfig_in->SetBranchAddress("Max", &StopMax);
  stopConfig_in->SetBranchAddress("Offset", &Offset);

  for (Int_t d_it = 0; d_it < NStops; ++d_it) {
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

  TH2 *MuonKinematics_seleff =
      SelMuExit ? GetHistogram<TH2>(EffCorrectorFile, "MuonKinematics_eff")
                : nullptr;
  TH2 *Ehadr_FV_detpos_seleff =
      GetHistogram<TH2>(EffCorrectorFile, "Ehadr_FV_detpos_seleff");

  EDep edr("EDeps", inpfile);

  TFile *of = CheckOpenFile(oupfile, "RECREATE");

  TH1D *EventRates_Selected =
      new TH1D("EventRates_Selected", ";Offset (cm);Count", 400, -250, 3750);
  TH1D *ERec_Selected = new TH1D(
      "ERec_Selected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count", 150, 0, 15);
  TH1D *EHadr_Selected =
      new TH1D("EHadr_Selected", ";E_{Hadr} (GeV);Count", 150, 0, 15);
  TH2D *EvRateERec_Selected = new TH2D(
      "EvRateERec_Selected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
      100, 0, 10, 400, -250, 3750);
  TH2D *EvRateEHadr_Selected =
      new TH2D("EvRateEHadr_Selected", ";E_{Hadr} (GeV);Offset (cm)", 100, 0,
               10, 400, -250, 3750);

  TH1D *EventRates_Corrected =
      new TH1D("EventRates_Corrected", ";Offset (cm);Count", 400, -250, 3750);
  TH1D *ERec_Corrected =
      new TH1D("ERec_Corrected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count",
               150, 0, 15);
  TH1D *EHadr_Corrected =
      new TH1D("EHadr_Corrected", ";E_{Hadr} (GeV);Count", 150, 0, 15);
  TH2D *EvRateERec_Corrected = new TH2D(
      "EvRateERec_Corrected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
      100, 0, 10, 400, -250, 3750);
  TH2D *EvRateEHadr_Corrected =
      new TH2D("EvRateEHadr_Corrected", ";E_{Hadr} (GeV);Offset (cm)", 100, 0,
               10, 400, -250, 3750);

  TDirectory *oupD = of;

  TDirectory *wD = oupD->mkdir("StopEventRates");
  wD->cd();

  std::vector<TH1D *> StopEventRates_Selected;
  std::vector<TH1D *> StopEventRates_Corrected;

  for (Int_t stop_it = 0; stop_it < NStops; ++stop_it) {
    StopEventRates_Selected.push_back(new TH1D(
        (std::string("StopEventRates_Selected_Stop") + to_str(stop_it)).c_str(),
        ";Offset (cm);Count", 400, -250, 3750));
    StopEventRates_Corrected.push_back(new TH1D(
        (std::string("StopEventRates_Corrected_Stop") + to_str(stop_it))
            .c_str(),
        ";Offset (cm);Count", 400, -250, 3750));
  }
  oupD->cd();

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

    if ((edr.stop < 0) || (SelMuExitKE && (edr.LepExitKE < SelMuExitKE)) ||
        (edr.TotalNonlep_Dep_veto > HadrVeto) || (edr.PrimaryLepPDG != 13)) {
      friendtree->Fill();
      continue;
    }

    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << edr.vtx[0] << ", " << edr.vtx[1] << ", " << edr.vtx[2]
                << "}, Enu: " << edr.nu_4mom[3] << " )" << std::endl;
    }

    if (SelMuExit) {
      TVector3 TrStr(edr.vtx[0], edr.vtx[1], edr.vtx[2]);
      TVector3 TrDir(edr.PrimaryLep_4mom[0], edr.PrimaryLep_4mom[1],
                     edr.PrimaryLep_4mom[2]);
      TrDir = TrDir.Unit();

      double MuToWall = CalculateToWall(BBs[edr.stop], TrStr, TrDir) * 1E-2;
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

    Int_t hadr_eff_x_bin = Ehadr_FV_detpos_seleff->GetXaxis()->FindFixBin(
        edr.ERecProxy_True - edr.PrimaryLep_4mom[3]);
    Int_t hadr_eff_y_bin =
        Ehadr_FV_detpos_seleff->GetYaxis()->FindFixBin(edr.vtxInDetX);

    if (Ehadr_FV_detpos_seleff->GetBinContent(hadr_eff_x_bin, hadr_eff_y_bin)) {
      hadreffweight =
          1.0 /
          Ehadr_FV_detpos_seleff->GetBinContent(hadr_eff_x_bin, hadr_eff_y_bin);
    } else {
      hadreffweight = 1;
    }

    effweight = (SelMuExit ? mueffweight : 1) * hadreffweight;

    // Fill selected histos
    EventRates_Selected->Fill(-1 * edr.vtx[0], edr.stop_weight);
    ERec_Selected->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                            edr.TotalNonlep_Dep_veto,
                        edr.stop_weight);
    EHadr_Selected->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                         edr.stop_weight);

    EvRateERec_Selected->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                                  edr.TotalNonlep_Dep_veto,
                              -1 * edr.vtx[0], edr.stop_weight);
    EvRateEHadr_Selected->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                               -1 * edr.vtx[0], edr.stop_weight);
    StopEventRates_Selected[edr.stop]->Fill(-1 * edr.vtx[0], edr.stop_weight);

    EventRates_Corrected->Fill(-1 * edr.vtx[0], edr.stop_weight * effweight);
    ERec_Corrected->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                             edr.TotalNonlep_Dep_veto,
                         edr.stop_weight * effweight);
    EHadr_Corrected->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                          edr.stop_weight * effweight);

    EvRateERec_Corrected->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                                   edr.TotalNonlep_Dep_veto,
                               -1 * edr.vtx[0], edr.stop_weight * effweight);
    EvRateEHadr_Corrected->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                                -1 * edr.vtx[0], edr.stop_weight * effweight);
    StopEventRates_Corrected[edr.stop]->Fill(-1 * edr.vtx[0],
                                             edr.stop_weight * effweight);

    friendtree->Fill();
  }

  of->Write();
  of->Close();
}
