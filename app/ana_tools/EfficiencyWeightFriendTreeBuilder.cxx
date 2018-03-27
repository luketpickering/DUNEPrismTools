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
bool CorrectInAbsPos = false;
std::vector<double> ERecBinning;

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
         "\t-A                          : Determine correction from absolute "
         "position and EProxy rather than hadronic shower kinematics.\n"
         "\t-b <binning descriptor>     : Energy binning descriptor that can "
         "be used to produce a perfect efficiency correction in the relevant "
         "off-axis bins.\n"
         "\t-FDproton                   : Build missing proton energ fake "
         "data. This is very hard coded, if you don't know what it is, don't "
         "use it.\n"
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
    } else if (std::string(argv[opt]) == "-E") {
      EffCorrectorFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-v") {
      HadrVeto = str2T<double>(argv[++opt]) * 1E-3;
    } else if (std::string(argv[opt]) == "-m") {
      SelMuExitKE = str2T<double>(argv[++opt]) * 1E-3;
      SelMuExit = true;
    } else if (std::string(argv[opt]) == "-A") {
      CorrectInAbsPos = true;
    } else if (std::string(argv[opt]) == "-b") {
      std::vector<std::string> binDescriptors =
          ParseToVect<std::string>(argv[++opt], ",");
      ERecBinning.clear();
      for (size_t vbd_it = 0; vbd_it < binDescriptors.size(); ++vbd_it) {
        AppendVect(ERecBinning, BuildDoubleList(binDescriptors[vbd_it]));
      }

      for (size_t bin_it = 1; bin_it < ERecBinning.size(); ++bin_it) {
        if (ERecBinning[bin_it] == ERecBinning[bin_it - 1]) {
          std::cout << "[INFO]: Removing duplciate bin low edge " << bin_it
                    << " low edge: " << ERecBinning[bin_it] << std::endl;
          ERecBinning.erase(ERecBinning.begin() + bin_it);
        }
      }

      for (size_t bin_it = 1; bin_it < ERecBinning.size(); ++bin_it) {
        if (ERecBinning[bin_it] < ERecBinning[bin_it - 1]) {
          std::cout << "[ERROR]: Bin " << bin_it
                    << " low edge: " << ERecBinning[bin_it]
                    << " is smaller than bin " << (bin_it - 1)
                    << " low edge: " << ERecBinning[bin_it - 1] << std::endl;
          exit(1);
        }
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

  if (!ERecBinning.size()) {
    int argc_dum = 3;
    char const *argv_dum[] = {"", "-b", "0_10:0.1"};
    handleOpts(argc_dum, argv_dum);
  }

  TChain *config_in = OpenTChainWithFileList("configTree", inpfile);
  Int_t NStops;
  Double_t FVGap[3];
  config_in->SetBranchAddress("NStops", &NStops);
  config_in->SetBranchAddress("FVGap", &FVGap);
  config_in->GetEntry(0);

  TChain *stopConfig_in = OpenTChainWithFileList("stopConfigTree", inpfile);
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
              << Min[1] << ", " << Min[2] << "} -- {" << Max[0] << ", "
              << Max[1] << ", " << Max[2] << "}" << std::endl;
  }

  TH2 *MuonKinematics_seleff =
      SelMuExit ? GetHistogram<TH2>(EffCorrectorFile, "MuonKinematics_eff")
                : nullptr;
  TH2 *Ehadr_FV_detpos_seleff =
      GetHistogram<TH2>(EffCorrectorFile, "Ehadr_FV_detpos_seleff");

  TH2 *Eproxy_FV_abspos_seleff =
      CorrectInAbsPos
          ? GetHistogram<TH2>(EffCorrectorFile, "Eproxy_FV_abspos_seleff")
          : nullptr;

  EDep edr("EDeps", inpfile);

  TFile *of = CheckOpenFile(oupfile, "RECREATE");

  TH1D *EventRates_Selected = nullptr;
  TH1D *ERec_Selected = nullptr;
  TH1D *EHadr_Selected = nullptr;
  TH2D *EvRateERec_Selected = nullptr;
  TH2D *EvRateEHadr_Selected = nullptr;
  TH2D *EvRateEProxy_Selected = nullptr;
  TH1D *EventRates_Corrected = nullptr;
  TH1D *ERec_Corrected = nullptr;
  TH1D *EHadr_Corrected = nullptr;
  TH2D *EvRateERec_Corrected = nullptr;
  TH2D *EvRateEHadr_Corrected = nullptr;
  TH2D *EvRateEProxy_Corrected = nullptr;

  std::vector<double> XRangeBins;
  if (CorrectInAbsPos) {
    for (Int_t bi_it = 0;
         bi_it < Eproxy_FV_abspos_seleff->GetYaxis()->GetNbins(); ++bi_it) {
      XRangeBins.push_back(
          Eproxy_FV_abspos_seleff->GetYaxis()->GetBinLowEdge(bi_it + 1));
    }
    XRangeBins.push_back(Eproxy_FV_abspos_seleff->GetYaxis()->GetBinUpEdge(
        Eproxy_FV_abspos_seleff->GetYaxis()->GetNbins()));

    EventRates_Selected = new TH1D("EventRates_Selected", ";Offset (cm);Count",
                                   (XRangeBins.size() - 1), XRangeBins.data());
    ERec_Selected =
        new TH1D("ERec_Selected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count",
                 (ERecBinning.size() - 1), ERecBinning.data());
    EHadr_Selected = new TH1D("EHadr_Selected", ";E_{Hadr} (GeV);Count",
                              (ERecBinning.size() - 1), ERecBinning.data());
    EvRateERec_Selected =
        new TH2D("EvRateERec_Selected",
                 ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
    EvRateEHadr_Selected =
        new TH2D("EvRateEHadr_Selected", ";E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
    EvRateEProxy_Selected =
        new TH2D("EvRateEProxy_Selected", ";E_{#nu, proxy} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
    EventRates_Corrected =
        new TH1D("EventRates_Corrected", ";Offset (cm);Count",
                 (XRangeBins.size() - 1), XRangeBins.data());
    ERec_Corrected =
        new TH1D("ERec_Corrected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count",
                 (ERecBinning.size() - 1), ERecBinning.data());
    EHadr_Corrected = new TH1D("EHadr_Corrected", ";E_{Hadr} (GeV);Count",
                               (ERecBinning.size() - 1), ERecBinning.data());
    EvRateERec_Corrected =
        new TH2D("EvRateERec_Corrected",
                 ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
    EvRateEHadr_Corrected =
        new TH2D("EvRateEHadr_Corrected", ";E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
    EvRateEProxy_Corrected =
        new TH2D("EvRateEProxy_Corrected", ";E_{#nu, proxy} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(),
                 (XRangeBins.size() - 1), XRangeBins.data());
  } else {
    EventRates_Selected =
        new TH1D("EventRates_Selected", ";Offset (cm);Count", 400, -250, 3750);
    ERec_Selected =
        new TH1D("ERec_Selected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count",
                 (ERecBinning.size() - 1), ERecBinning.data());
    EHadr_Selected = new TH1D("EHadr_Selected", ";E_{Hadr} (GeV);Count",
                              (ERecBinning.size() - 1), ERecBinning.data());
    EvRateERec_Selected =
        new TH2D("EvRateERec_Selected",
                 ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(), 400, -250, 3750);
    EvRateEHadr_Selected =
        new TH2D("EvRateEHadr_Selected", ";E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(), 400, -250, 3750);
    EvRateEProxy_Selected =
        new TH2D("EvRateEProxy_Selected", ";E_{#nu, proxy} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(), 400, -250, 3750);
    EventRates_Corrected =
        new TH1D("EventRates_Corrected", ";Offset (cm);Count", 400, -250, 3750);
    ERec_Corrected =
        new TH1D("ERec_Corrected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count",
                 (ERecBinning.size() - 1), ERecBinning.data());
    EHadr_Corrected = new TH1D("EHadr_Corrected", ";E_{Hadr} (GeV);Count",
                               (ERecBinning.size() - 1), ERecBinning.data());
    EvRateERec_Corrected =
        new TH2D("EvRateERec_Corrected",
                 ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(), 400, -250, 3750);
    EvRateEHadr_Corrected =
        new TH2D("EvRateEHadr_Corrected", ";E_{Hadr} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(), 400, -250, 3750);
    EvRateEProxy_Corrected =
        new TH2D("EvRateEProxy_Corrected", ";E_{#nu, proxy} (GeV);Offset (cm)",
                 (ERecBinning.size() - 1), ERecBinning.data(), 400, -250, 3750);
  }

  TDirectory *oupD = of;

  TDirectory *sliceD = oupD->mkdir("SliceEventRates");
  sliceD->cd();

  std::vector<TH1D *> SliceERec_Selected;
  std::vector<TH1D *> SliceERec_Corrected;
  TH1D *SliceHelper = nullptr;
  if (CorrectInAbsPos) {
    SliceHelper =
        new TH1D("SliceHelper", "", (XRangeBins.size() - 1), XRangeBins.data());
    SliceHelper->SetDirectory(nullptr);

    for (size_t bin_it = 1; bin_it < XRangeBins.size(); ++bin_it) {
      SliceERec_Selected.push_back(new TH1D(
          (std::string("SliceERec_Selected") + to_str(bin_it - 1)).c_str(),
          ";E_{Rec} (GeV);Count", 100, 0, 10));
      SliceERec_Corrected.push_back(new TH1D(
          (std::string("SliceERec_Corrected") + to_str(bin_it - 1)).c_str(),
          ";E_{Rec} (GeV);Count", 100, 0, 10));
    }
  }

  TDirectory *stopD = oupD->mkdir("StopEventRates");
  stopD->cd();

  std::vector<TH1D *> StopEventRates_Selected;
  std::vector<TH1D *> StopEventRates_Corrected;

  std::vector<TH1D *> StopEProxy_Selected;
  std::vector<TH1D *> StopEProxy_Corrected;

  std::vector<TH1D *> StopERec_Selected;
  std::vector<TH1D *> StopERec_Corrected;

  for (Int_t stop_it = 0; stop_it < NStops; ++stop_it) {
    StopEventRates_Selected.push_back(new TH1D(
        (std::string("StopEventRates_Selected_Stop") + to_str(stop_it)).c_str(),
        ";Offset (cm);Count", 400, -250, 3750));
    StopEventRates_Corrected.push_back(new TH1D(
        (std::string("StopEventRates_Corrected_Stop") + to_str(stop_it))
            .c_str(),
        ";Offset (cm);Count", 400, -250, 3750));

    StopEProxy_Selected.push_back(new TH1D(
        (std::string("StopEProxy_Selected_Stop") + to_str(stop_it)).c_str(),
        ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));
    StopEProxy_Corrected.push_back(new TH1D(
        (std::string("StopEProxy_Corrected_Stop") + to_str(stop_it)).c_str(),
        ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));

    StopERec_Selected.push_back(new TH1D(
        (std::string("StopERec_Selected_Stop") + to_str(stop_it)).c_str(),
        ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));
    StopERec_Corrected.push_back(new TH1D(
        (std::string("StopERec_Corrected_Stop") + to_str(stop_it)).c_str(),
        ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));
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

    // if (BuildMissingProtonEFakeData) {
    //   ProtonFakeDataWeight_outputtree->GetEntry(e_it);

    //   edr.stop_weight *= ProtonFakeDataWeight;
    //   edr.TotalNonlep_Dep_veto -= edr.ProtonDep_veto * 0.2;
    //   edr.TotalNonlep_Dep_FV -= edr.ProtonDep_FV * 0.2;
    //   edr.ERecProxy_True -= edr.EKinProton_True * 0.2;
    // }

    if ((edr.stop < 0) || (SelMuExitKE && (edr.LepExitKE < SelMuExitKE)) ||
        (edr.TotalNonlep_Dep_veto > HadrVeto) || (edr.PrimaryLepPDG != 13)) {
      std::cout << "[ERROR]: Currently EfficiencyWeightFriendTreeBuilder "
                   "expected to receive pre-selected events."
                << std::endl;
      exit(1);
      friendtree->Fill();
      continue;
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

    if (!CorrectInAbsPos) {
      Int_t hadr_eff_x_bin = Ehadr_FV_detpos_seleff->GetXaxis()->FindFixBin(
          edr.ERecProxy_True - edr.PrimaryLep_4mom[3]);
      Int_t hadr_eff_y_bin =
          Ehadr_FV_detpos_seleff->GetYaxis()->FindFixBin(edr.vtxInDetX);

      if (Ehadr_FV_detpos_seleff->GetBinContent(hadr_eff_x_bin,
                                                hadr_eff_y_bin)) {
        hadreffweight = 1.0 /
                        Ehadr_FV_detpos_seleff->GetBinContent(hadr_eff_x_bin,
                                                              hadr_eff_y_bin);
      } else {
        hadreffweight = 1;
      }
      effweight = hadreffweight;
    } else {
      Int_t hadr_eff_x_bin =
          Eproxy_FV_abspos_seleff->GetXaxis()->FindFixBin(edr.ERecProxy_True);
      Int_t hadr_eff_y_bin =
          Eproxy_FV_abspos_seleff->GetYaxis()->FindFixBin(-1 * edr.vtx[0]);

      if (Eproxy_FV_abspos_seleff->GetBinContent(hadr_eff_x_bin,
                                                 hadr_eff_y_bin)) {
        hadreffweight = 1.0 /
                        Eproxy_FV_abspos_seleff->GetBinContent(hadr_eff_x_bin,
                                                               hadr_eff_y_bin);
      } else {
        hadreffweight = 1;
      }
      effweight = (SelMuExit ? mueffweight : 1) * hadreffweight;

      Int_t xb = SliceHelper->GetXaxis()->FindFixBin(-1 * edr.vtx[0]);
      if (xb && (xb <= Int_t(SliceERec_Selected.size()))) {
        SliceERec_Selected[xb - 1]->Fill(edr.PrimaryLep_4mom[3] +
                                             edr.TotalNonlep_Dep_FV +
                                             edr.TotalNonlep_Dep_veto,
                                         edr.stop_weight);
        SliceERec_Corrected[xb - 1]->Fill(edr.PrimaryLep_4mom[3] +
                                              edr.TotalNonlep_Dep_FV +
                                              edr.TotalNonlep_Dep_veto,
                                          edr.stop_weight * effweight);
      }
    }

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
    EvRateEProxy_Selected->Fill(edr.ERecProxy_True, -1 * edr.vtx[0],
                                edr.stop_weight);
    StopEProxy_Selected[edr.stop]->Fill(edr.ERecProxy_True, edr.stop_weight);
    StopERec_Selected[edr.stop]->Fill(edr.PrimaryLep_4mom[3] +
                                          edr.TotalNonlep_Dep_FV +
                                          edr.TotalNonlep_Dep_veto,
                                      edr.stop_weight);

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
    EvRateEProxy_Corrected->Fill(edr.ERecProxy_True, -1 * edr.vtx[0],
                                 edr.stop_weight * effweight);
    StopEProxy_Corrected[edr.stop]->Fill(edr.ERecProxy_True,
                                         edr.stop_weight * effweight);
    StopERec_Corrected[edr.stop]->Fill(edr.PrimaryLep_4mom[3] +
                                           edr.TotalNonlep_Dep_FV +
                                           edr.TotalNonlep_Dep_veto,
                                       edr.stop_weight * effweight);
    friendtree->Fill();
  }

  of->Write();
  of->Close();
}
