#include "EDepTreeReader.h"
#include "TTree.h"
#include "Utils.hxx"

#include <string>
#include <vector>

std::string inpfile;
std::string oupfile;
std::string lincombfile;
bool BuildMissingProtonEFakeData = false;

double HadrVeto = 0;  // GeV
bool SelMuExit = false;
double SelMuExitKE = 0;  // GeV
// std::vector<double> ExtraFV;

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "\t-i <fulldetprocess.root>    : TChain descriptor for"
         " input tree. \n"
         "\t-F <fluxfitresults.root>             : Result of "
         "dp_FitFluxes to use to build observation.\n"
         "\t-o <outputfile.root>        : Output file to write "
         "selected tree to.\n"
         "\t-v <hadr veto threshold>    : Hadronic shower veto threshold in "
         "MeV.\n"
         "\t-FDproton                   : Build missing proton energy fake "
         "data. This is very hard coded, if you don't know what it is, don't "
         "use it.\n"
         // "\t-XFV <fvx,y,z>              : Extra FV reduction applied at the "
         // "selection level.\n"
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
    } else if (std::string(argv[opt]) == "-F") {
      lincombfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      oupfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-v") {
      HadrVeto = str2T<double>(argv[++opt]) * 1E-3;
    } else if (std::string(argv[opt]) == "-m") {
      SelMuExitKE = str2T<double>(argv[++opt]) * 1E-3;
      SelMuExit = true;
    } else if (std::string(argv[opt]) == "-FDproton") {
      BuildMissingProtonEFakeData = true;
    }
    // else if (std::string(argv[opt]) == "-XFV") {
    //   ExtraFV = ParseToVect<double>(argv[++opt], ",");
    // }
    else {
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

  // if(!ExtraFV.size()){
  //   ExtraFV.push_back(0);
  //   ExtraFV.push_back(0);
  //   ExtraFV.push_back(0);
  // }

  TChain *config_in = OpenTChainWithFileList("configTree", inpfile);
  Int_t NStops;
  config_in->SetBranchAddress("NStops", &NStops);
  double timesep_us = 0xdeadbeef;
  if (CheckTTreeHasBranch(config_in, "timesep_us")) {
    config_in->SetBranchAddress("timesep_us", &timesep_us);
  }
  config_in->GetEntry(0);

  TH1D *CoeffWeightingHelper = nullptr;
  TH2D *OffSetVsEProxy = nullptr;
  TH2D *OffSetVsETrue = nullptr;
  std::vector<TH1D *> SliceEProxy_True_LinCombWeighted;
  std::vector<TH1D *> SliceETrue_True;
  std::vector<TH1D *> SliceERec_Selected;
  std::vector<TH1D *> SliceEMu_Selected;
  std::vector<TH1D *> SliceEMu_True;
  std::vector<TH1D *> SliceEHadr_Selected;
  std::vector<TH1D *> SliceETHadr_Selected;

  if (lincombfile.size()) {
    TChain *CoeffTree = new TChain("CoeffTree");

    CoeffTree->Add(lincombfile.c_str());

    double XRange[2];
    double Coeff;

    CoeffTree->SetBranchAddress("XRange", &XRange);
    CoeffTree->SetBranchAddress("Coeff", &Coeff);

    std::vector<double> XRangeBins;
    std::vector<double> Coeffs;

    std::cout << "[INFO]: XRange bins: " << std::flush;
    CoeffTree->GetEntry(0);
    XRangeBins.push_back(XRange[0]);
    Coeffs.push_back(Coeff);
    XRangeBins.push_back(XRange[1]);
    std::cout << XRangeBins[0] << ", " << XRangeBins[1] << std::flush;
    for (Long64_t i = 1; i < CoeffTree->GetEntries(); ++i) {
      CoeffTree->GetEntry(i);

      // If non-contiguous, must push an empty bit between.
      if (fabs(XRangeBins.back() - XRange[0]) > 1E-5) {
        Coeffs.push_back(0);
        XRangeBins.push_back(XRange[0]);
        std::cout << ", " << XRangeBins.back() << std::flush;
      }

      Coeffs.push_back(Coeff);
      XRangeBins.push_back(XRange[1]);
      std::cout << ", " << XRangeBins.back() << std::flush;
    }
    std::cout << std::endl;

    CoeffWeightingHelper = new TH1D("CoeffWeightingHelper", "",
                                    (XRangeBins.size() - 1), XRangeBins.data());
    CoeffWeightingHelper->SetDirectory(nullptr);

    for (size_t bin_it = 1; bin_it < XRangeBins.size(); ++bin_it) {
      CoeffWeightingHelper->SetBinContent(bin_it, Coeffs[bin_it - 1]);

      SliceEProxy_True_LinCombWeighted.push_back(
          new TH1D((std::string("SliceEProxy_True_LinCombWeighted_Slice") +
                    to_str(bin_it - 1))
                       .c_str(),
                   ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));

      SliceETrue_True.push_back(new TH1D(
          (std::string("SliceETrue_True") + to_str(bin_it - 1)).c_str(),
          ";E_{#nu} (GeV);Count", 100, 0, 10));
      SliceERec_Selected.push_back(new TH1D(
          (std::string("SliceERec_Selected") + to_str(bin_it - 1)).c_str(),
          ";E_{Rec} (GeV);Count", 100, 0, 10));
      SliceEMu_Selected.push_back(new TH1D(
          (std::string("SliceEMu_Selected") + to_str(bin_it - 1)).c_str(),
          ";E_{#mu} (GeV);Count", 100, 0, 10));
      SliceEMu_True.push_back(
          new TH1D((std::string("SliceEMu_True") + to_str(bin_it - 1)).c_str(),
                   ";E_{#mu} (GeV);Count", 100, 0, 10));
      SliceEHadr_Selected.push_back(new TH1D(
          (std::string("SliceEHadr_Selected") + to_str(bin_it - 1)).c_str(),
          ";E_{Hadr} (GeV);Count", 100, 0, 10));
      SliceETHadr_Selected.push_back(new TH1D(
          (std::string("SliceETHadr_Selected") + to_str(bin_it - 1)).c_str(),
          ";E_{Hadr} (GeV);Count", 100, 0, 10));
    }

    OffSetVsEProxy =
        new TH2D("OffSetVsEProxy", ";E_{#nu,Proxy} (GeV);Offset (cm);Count",
                 100, 0, 10, XRangeBins.size() - 1, XRangeBins.data());
    OffSetVsETrue =
        new TH2D("OffSetVsETrue", ";E_{#nu} (GeV);Offset (cm);Count", 100, 0,
                 10, XRangeBins.size() - 1, XRangeBins.data());

    OffSetVsEProxy->SetDirectory(nullptr);
    OffSetVsETrue->SetDirectory(nullptr);
  }

  TChain *stopConfig_in = OpenTChainWithFileList("stopConfigTree", inpfile);

  EDep edr("EDeps", inpfile);

  Double_t ProtonFakeDataWeight = 1;
  if (BuildMissingProtonEFakeData) {
    // Should probably check that this exists, but leaving hardcoded and
    // presumptive for now as CheckTTreeHasBranch doesn't seem to work on
    // pre-friended trees.
    edr.tree->SetBranchAddress("EnuTp", &ProtonFakeDataWeight);
  }

  TFile *of = CheckOpenFile(oupfile, "RECREATE");
  of->cd();

  TTree *config_copy = config_in->CloneTree();
  TTree *stopConfig_copy = stopConfig_in->CloneTree();
  config_copy->SetDirectory(of);
  stopConfig_copy->SetDirectory(of);

  TTree *selconftree = new TTree("SelectionConfigTree", "");
  Int_t NTotal = 0, NInStops = 0, NNumuCC = 0, NNumuNC = 0, NNumubCC = 0,
        NNumubNC = 0, NNueCC = 0, NNueNC = 0, NNuebCC = 0, NNuebNC = 0,
        NSel = 0;
  Bool_t SelectOnMuonExit = SelMuExit;
  Double_t MuonExitKECut_MeV = SelMuExitKE * 1E3,
           HadronicShowerVetoCut_MeV = HadrVeto * 1E3;
  selconftree->Branch("NTotal", &NTotal, "NTotal/I");
  selconftree->Branch("NInStops", &NInStops, "NInStops/I");
  selconftree->Branch("NNumuCC", &NNumuCC, "NNumuCC/I");
  selconftree->Branch("NNumuNC", &NNumuNC, "NNumuNC/I");
  selconftree->Branch("NNumubCC", &NNumubCC, "NNumubCC/I");
  selconftree->Branch("NNumubNC", &NNumubNC, "NNumubNC/I");
  selconftree->Branch("NNueCC", &NNueCC, "NNueCC/I");
  selconftree->Branch("NNueNC", &NNueNC, "NNueNC/I");
  selconftree->Branch("NNuebCC", &NNuebCC, "NNuebCC/I");
  selconftree->Branch("NNuebNC", &NNuebNC, "NNuebNC/I");
  selconftree->Branch("NSel", &NSel, "NSel/I");
  selconftree->Branch("SelectOnMuonExit", &SelectOnMuonExit,
                      "SelectOnMuonExit/O");
  selconftree->Branch("MuonExitKECut_MeV", &MuonExitKECut_MeV,
                      "MuonExitKECut_MeV/D");
  selconftree->Branch("HadronicShowerVetoCut_MeV", &HadronicShowerVetoCut_MeV,
                      "HadronicShowerVetoCut_MeV/D");

  TTree *OutputTree = new TTree("EDeps", "");

  TH1D *EventRates_Selected =
      new TH1D("EventRates_Selected", ";Offset (cm);Count", 400, -250, 3750);
  TH1D *EventRates_True =
      new TH1D("EventRates_True", ";Offset (cm);Count", 400, -250, 3750);

  TH1D *ERec_Selected = new TH1D(
      "ERec_Selected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Count", 100, 0, 10);
  TH1D *ETrue_Selected =
      new TH1D("ETrue_Selected", ";E_{#nu} (GeV);Count", 100, 0, 10);
  TH1D *EHadr_Selected =
      new TH1D("EHadr_Selected", ";E_{Hadr} (GeV);Count", 100, 0, 10);
  TH1D *EHadr_True =
      new TH1D("EHadr_True", ";E_{Hadr} (GeV);Count", 100, 0, 10);

  TH2D *EvRateERec_Selected = new TH2D(
      "EvRateERec_Selected", ";E_{Rec} = E_{#mu} + E_{Hadr} (GeV);Offset (cm)",
      100, 0, 10, 400, -250, 3750);

  TH2D *EvRateEHadr_Selected =
      new TH2D("EvRateEHadr_Selected", ";E_{Hadr} (GeV);Offset (cm)", 100, 0,
               10, 400, -250, 3750);
  TH2D *EvRateEHadr_True =
      new TH2D("EvRateEHadr_True", ";E_{Hadr} (GeV);Offset (cm)", 100, 0, 10,
               400, -250, 3750);

  TH2D *EvRateEProxy_Selected =
      new TH2D("EvRateEProxy_Selected", ";E_{#nu,proxy} (GeV);Offset (cm)", 100,
               0, 10, 400, -250, 3750);
  TH2D *EvRateEProxy_True =
      new TH2D("EvRateEProxy_True", ";E_{#nu,proxy} (GeV);Offset (cm)", 100, 0,
               10, 400, -250, 3750);

  TH2D *EvRateEVeto = new TH2D("EvRateEVeto", ";E_{Hadr} (GeV);Offset (cm)",
                               100, 0, 2, 400, -250, 3750);

  TH2D *EProxyERec_Selected =
      new TH2D("EProxyERec_Selected", ";E_{#nu,proxy} (GeV);E_{Rec} (GeV)", 100,
               1E-5, 10, 100, 1E-5, 10);
  TH2D *ETrueEProxy_Selected =
      new TH2D("ETrueEProxy_Selected", ";E_{#nu} (GeV);E_{#nu,proxy} (GeV)",
               100, 1E-5, 10, 100, 1E-5, 10);
  TH2D *ETrueERec_Selected =
      new TH2D("ETrueERec_Selected", ";E_{#nu} (GeV);E_{Rec} (GeV)", 100, 1E-5,
               10, 100, 1E-5, 10);

  TH2D *ETrueEmu_Selected =
      new TH2D("ETrueEmu_Selected", ";E_{#nu} (GeV);E_{#mu} (GeV)", 100, 1E-5,
               10, 100, 1E-5, 10);

  TH2D *ETrueERecNeutron_Selected = new TH2D(
      "ETrueERecNeutron_Selected", ";E_{#nu} (GeV);E_{Rec} Neutron (GeV)", 100,
      1E-5, 10, 100, 1E-5, 10);
  TH2D *ETrueERecNeutron_timesep_Selected = new TH2D(
      "ETrueERecNeutron_timesep_Selected",
      ";E_{#nu} (GeV);E_{Rec} Neutron (GeV)", 100, 1E-5, 10, 100, 1E-5, 10);

  TH2D *ETrueERecProton_Selected = new TH2D(
      "ETrueERecProton_Selected", ";E_{#nu} (GeV);E_{Rec} Proton (GeV)", 100,
      1E-5, 10, 100, 1E-5, 10);
  TH2D *ETrueERecProton_timesep_Selected = new TH2D(
      "ETrueERecProton_timesep_Selected", ";E_{#nu} (GeV);E_{Rec} Proton (GeV)",
      100, 1E-5, 10, 100, 1E-5, 10);

  TH2D *ETrueERecPion_Selected =
      new TH2D("ETrueERecPion_Selected", ";E_{#nu} (GeV);E_{Rec} Pion (GeV)",
               100, 1E-5, 10, 100, 1E-5, 10);
  TH2D *ETrueERecPion_timesep_Selected = new TH2D(
      "ETrueERecPion_timesep_Selected", ";E_{#nu} (GeV);E_{Rec} Pion (GeV)",
      100, 1E-5, 10, 100, 1E-5, 10);

  TH2D *ETrueERecOther_Selected =
      new TH2D("ETrueERecOther_Selected", ";E_{#nu} (GeV);E_{Rec} Other (GeV)",
               100, 1E-5, 10, 100, 1E-5, 10);
  TH2D *ETrueERecOther_timesep_Selected = new TH2D(
      "ETrueERecOther_timesep_Selected", ";E_{#nu} (GeV);E_{Rec} Other (GeV)",
      100, 1E-5, 10, 100, 1E-5, 10);

  TDirectory *oupD = of;

  TDirectory *wD = oupD->mkdir("StopEventRates");
  wD->cd();

  std::vector<TH1D *> StopEventRates_Selected;
  std::vector<TH1D *> StopEProxy_Selected;
  std::vector<TH1D *> StopEProxy_True;
  std::vector<TH1D *> StopETrue_Selected;
  std::vector<TH1D *> StopETrue_True;
  std::vector<TH1D *> StopERec_Selected;
  std::vector<TH1D *> StopEventRates_True;
  std::vector<TH1D *> StopEventRates_TrueUnWeighted;

  for (Int_t stop_it = 0; stop_it < NStops; ++stop_it) {
    StopEventRates_Selected.push_back(new TH1D(
        (std::string("StopEventRates_Selected_Stop") + to_str(stop_it)).c_str(),
        ";Offset (cm);Count", 400, -250, 3750));
    StopEProxy_Selected.push_back(new TH1D(
        (std::string("StopEProxy_Selected_Stop") + to_str(stop_it)).c_str(),
        ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));
    StopEProxy_True.push_back(new TH1D(
        (std::string("StopEProxy_True_Stop") + to_str(stop_it)).c_str(),
        ";E_{#nu, proxy} (GeV);Count", 100, 0, 10));
    StopEventRates_True.push_back(new TH1D(
        (std::string("StopEventRates_True_Stop") + to_str(stop_it)).c_str(),
        ";Offset (cm);Count", 400, -250, 3750));
    StopEventRates_TrueUnWeighted.push_back(new TH1D(
        (std::string("StopEventRates_TrueUnWeighted_Stop") + to_str(stop_it))
            .c_str(),
        ";Offset (cm);Count", 400, -250, 3750));

    StopETrue_Selected.push_back(
        new TH1D((std::string("StopETrue_Selected") + to_str(stop_it)).c_str(),
                 ";E_{#nu} (GeV);Count", 100, 0, 10));
    StopETrue_True.push_back(
        new TH1D((std::string("StopETrue_True") + to_str(stop_it)).c_str(),
                 ";E_{#nu} (GeV);Count", 100, 0, 10));
    StopERec_Selected.push_back(
        new TH1D((std::string("StopERec_Selected") + to_str(stop_it)).c_str(),
                 ";E_{Rec} (GeV);Count", 100, 0, 10));
  }
  oupD->cd();

  EDep *OutputEDep = EDep::MakeTreeWriter(OutputTree, timesep_us, true);

  size_t loud_every = edr.GetEntries() / 10;
  Long64_t NEntries = edr.GetEntries();
  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    edr.GetEntry(e_it);

    NTotal++;
    NInStops += (edr.stop > -1);
    NNumuCC += (edr.stop > -1) && (edr.PrimaryLepPDG == 13);
    NNumuNC += (edr.stop > -1) && (edr.PrimaryLepPDG == 14);
    NNumubCC += (edr.stop > -1) && (edr.PrimaryLepPDG == -13);
    NNumubNC += (edr.stop > -1) && (edr.PrimaryLepPDG == -14);
    NNueCC += (edr.stop > -1) && (edr.PrimaryLepPDG == 11);
    NNueNC += (edr.stop > -1) && (edr.PrimaryLepPDG == 12);
    NNuebCC += (edr.stop > -1) && (edr.PrimaryLepPDG == -11);
    NNuebNC += (edr.stop > -1) && (edr.PrimaryLepPDG == -12);

    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << edr.vtx[0] << ", " << edr.vtx[1] << ", " << edr.vtx[2]
                << "}, Enu: " << edr.nu_4mom[3] << " )" << std::endl;
      if (BuildMissingProtonEFakeData) {
        std::cout << "\tProton fake data weight = " << ProtonFakeDataWeight
                  << std::endl;
      }
    }

    if ((edr.stop < 0) || (edr.PrimaryLepPDG != 13)) {
      continue;
    }

    if (BuildMissingProtonEFakeData) {
      edr.stop_weight *= ProtonFakeDataWeight;
      edr.TotalNonlep_Dep_veto -= edr.ProtonDep_veto * 0.2;
      edr.TotalNonlep_Dep_FV -= edr.ProtonDep_FV * 0.2;
      edr.ERecProxy_True -= edr.EKinProton_True * 0.2;
    }

    // Fill truth histos
    EventRates_True->Fill(-1 * edr.vtx[0], edr.stop_weight);
    EHadr_True->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                     edr.stop_weight);
    EvRateEHadr_True->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                           -1 * edr.vtx[0], edr.stop_weight);
    EvRateEProxy_True->Fill(edr.ERecProxy_True, -1 * edr.vtx[0],
                            edr.stop_weight);
    StopEProxy_True[edr.stop]->Fill(edr.ERecProxy_True, edr.stop_weight);
    StopEventRates_True[edr.stop]->Fill(-1 * edr.vtx[0], edr.stop_weight);
    StopEventRates_TrueUnWeighted[edr.stop]->Fill(-1 * edr.vtx[0]);
    EvRateEVeto->Fill(edr.TotalNonlep_Dep_veto, -1 * edr.vtx[0]);
    StopETrue_True[edr.stop]->Fill(edr.nu_4mom[3], edr.stop_weight);

    Int_t xb = -1;
    if (CoeffWeightingHelper) {
      OffSetVsEProxy->Fill(edr.ERecProxy_True, -1 * edr.vtx[0]);
      OffSetVsETrue->Fill(edr.nu_4mom[3], -1 * edr.vtx[0]);

      xb = CoeffWeightingHelper->GetXaxis()->FindFixBin(-1 * edr.vtx[0]);
      double CoeffWeight = CoeffWeightingHelper->GetBinContent(xb);

      if (xb && (xb <= Int_t(SliceEProxy_True_LinCombWeighted.size()))) {
        SliceEProxy_True_LinCombWeighted[xb - 1]->Fill(
            edr.ERecProxy_True, edr.stop_weight * CoeffWeight);

        SliceETrue_True[xb - 1]->Fill(edr.nu_4mom[3], edr.stop_weight);

        SliceEMu_True[xb - 1]->Fill(edr.PrimaryLep_4mom[3], edr.stop_weight);
      }
    }

    if ((edr.stop < 0) || (SelMuExitKE && (edr.LepExitKE < SelMuExitKE)) ||
        (edr.TotalNonlep_Dep_veto > HadrVeto) || (edr.PrimaryLepPDG != 13)) {
      continue;
    }

    NSel++;
    // Fill selected histos
    EventRates_Selected->Fill(-1 * edr.vtx[0], edr.stop_weight);
    ERec_Selected->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                            edr.TotalNonlep_Dep_veto,
                        edr.stop_weight);
    EHadr_Selected->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                         edr.stop_weight);
    ETrue_Selected->Fill(edr.nu_4mom[3], edr.stop_weight);

    EvRateERec_Selected->Fill(edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                                  edr.TotalNonlep_Dep_veto,
                              -1 * edr.vtx[0], edr.stop_weight);
    EvRateEHadr_Selected->Fill(edr.ERecProxy_True - edr.PrimaryLep_4mom[3],
                               -1 * edr.vtx[0], edr.stop_weight);
    EvRateEProxy_Selected->Fill(edr.ERecProxy_True, -1 * edr.vtx[0],
                                edr.stop_weight);
    StopEProxy_Selected[edr.stop]->Fill(edr.ERecProxy_True, edr.stop_weight);
    StopEventRates_Selected[edr.stop]->Fill(-1 * edr.vtx[0], edr.stop_weight);
    StopETrue_Selected[edr.stop]->Fill(edr.nu_4mom[3], edr.stop_weight);
    StopERec_Selected[edr.stop]->Fill(edr.PrimaryLep_4mom[3] +
                                          edr.TotalNonlep_Dep_FV +
                                          edr.TotalNonlep_Dep_veto,
                                      edr.stop_weight);

    EProxyERec_Selected->Fill(edr.ERecProxy_True,
                              edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                                  edr.TotalNonlep_Dep_veto,
                              edr.stop_weight);
    ETrueEProxy_Selected->Fill(edr.nu_4mom[3], edr.ERecProxy_True,
                               edr.stop_weight);
    ETrueERec_Selected->Fill(edr.nu_4mom[3],
                             edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                                 edr.TotalNonlep_Dep_veto,
                             edr.stop_weight);

    ETrueEmu_Selected->Fill(edr.nu_4mom[3], edr.PrimaryLep_4mom[3],
                            edr.stop_weight);
    ETrueERecNeutron_Selected->Fill(edr.nu_4mom[3],
                                    edr.NeutronDep_FV + edr.NeutronDep_veto,
                                    edr.stop_weight);
    ETrueERecNeutron_timesep_Selected->Fill(
        edr.nu_4mom[3], edr.NeutronDep_timesep_FV + edr.NeutronDep_timesep_veto,
        edr.stop_weight);
    ETrueERecProton_Selected->Fill(
        edr.nu_4mom[3], edr.ProtonDep_FV + edr.ProtonDep_veto, edr.stop_weight);
    ETrueERecProton_timesep_Selected->Fill(
        edr.nu_4mom[3], edr.ProtonDep_timesep_FV + edr.ProtonDep_timesep_veto,
        edr.stop_weight);
    ETrueERecPion_Selected->Fill(
        edr.nu_4mom[3],
        edr.PiCDep_FV + edr.PiCDep_veto + edr.Pi0Dep_FV + edr.Pi0Dep_veto,
        edr.stop_weight);
    ETrueERecPion_timesep_Selected->Fill(
        edr.nu_4mom[3], edr.PiCDep_timesep_FV + edr.PiCDep_timesep_veto +
                            edr.Pi0Dep_timesep_FV + edr.Pi0Dep_timesep_veto,
        edr.stop_weight);
    ETrueERecOther_Selected->Fill(
        edr.nu_4mom[3], edr.OtherDep_FV + edr.OtherDep_veto, edr.stop_weight);
    ETrueERecOther_timesep_Selected->Fill(
        edr.nu_4mom[3], edr.OtherDep_timesep_FV + edr.OtherDep_timesep_veto,
        edr.stop_weight);

    if (CoeffWeightingHelper) {
      if (xb && (xb <= Int_t(SliceEProxy_True_LinCombWeighted.size()))) {
        SliceERec_Selected[xb - 1]->Fill(edr.PrimaryLep_4mom[3] +
                                             edr.TotalNonlep_Dep_FV +
                                             edr.TotalNonlep_Dep_veto,
                                         edr.stop_weight);
        SliceEMu_Selected[xb - 1]->Fill(edr.PrimaryLep_4mom[3],
                                        edr.stop_weight);
        SliceEHadr_Selected[xb - 1]->Fill(
            edr.TotalNonlep_Dep_FV + edr.TotalNonlep_Dep_veto, edr.stop_weight);
        SliceETHadr_Selected[xb - 1]->Fill(
            edr.ERecProxy_True - edr.PrimaryLep_4mom[3], edr.stop_weight);
      }
    }

    OutputEDep->Copy(edr);
    OutputTree->Fill();
  }

  if (CoeffWeightingHelper) {
    TDirectory *oupD = of;

    TDirectory *wD = oupD->mkdir("TruthLinearCombinations");
    wD->cd();

    TH1D *LinCombEProxy = new TH1D(
        "LinCombEProxy", ";E_{#nu,Proxy} (GeV);Events / GeV", 100, 0, 10);
    TH1D *LinCombETrue =
        new TH1D("LinCombETrue", ";E_{#nu} (GeV);Events / GeV", 100, 0, 10);

    for (Int_t xbin = 1; xbin < OffSetVsETrue->GetXaxis()->GetNbins() + 1;
         ++xbin) {
      double bc_etrue = 0, be_etrue = 0, bc_eprox = 0, be_eprox = 0;

      for (Int_t ybin = 1; ybin < OffSetVsETrue->GetYaxis()->GetNbins() + 1;
           ++ybin) {
        bc_etrue += OffSetVsETrue->GetBinContent(xbin, ybin) *
                    CoeffWeightingHelper->GetBinContent(ybin);

        be_etrue += OffSetVsETrue->GetBinError(xbin, ybin) *
                    CoeffWeightingHelper->GetBinContent(ybin);

        bc_eprox += OffSetVsEProxy->GetBinContent(xbin, ybin) *
                    CoeffWeightingHelper->GetBinContent(ybin);
        be_eprox += OffSetVsEProxy->GetBinError(xbin, ybin) *
                    CoeffWeightingHelper->GetBinContent(ybin);
      }

      LinCombETrue->SetBinContent(xbin, bc_etrue);
      LinCombETrue->SetBinError(xbin, be_etrue);

      LinCombEProxy->SetBinContent(xbin, bc_eprox);
      LinCombEProxy->SetBinError(xbin, be_eprox);
    }

    CoeffWeightingHelper->SetDirectory(wD);
    OffSetVsEProxy->SetDirectory(wD);
    OffSetVsETrue->SetDirectory(wD);

    TDirectory *sliceD = oupD->mkdir("SliceObservables");

    for (size_t i = 0; i < SliceEProxy_True_LinCombWeighted.size(); ++i) {
      SliceEProxy_True_LinCombWeighted[i]->SetDirectory(wD);
      SliceETrue_True[i]->SetDirectory(sliceD);
      SliceERec_Selected[i]->SetDirectory(sliceD);
      SliceEMu_True[i]->SetDirectory(sliceD);
      SliceEMu_Selected[i]->SetDirectory(sliceD);
      SliceEHadr_Selected[i]->SetDirectory(sliceD);
      SliceETHadr_Selected[i]->SetDirectory(sliceD);
    }
  }

  TH2D *EProxyERec_Selected_ETrueNorm = static_cast<TH2D *>(
      EProxyERec_Selected->Clone("EProxyERec_SelectedETrueNorm"));
  TH2D *ETrueEProxy_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueEProxy_Selected->Clone("ETrueEProxy_SelectedETrueNorm"));
  TH2D *ETrueERec_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueERec_Selected->Clone("ETrueERec_SelectedETrueNorm"));
  TH2D *ETrueEmu_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueEmu_Selected->Clone("ETrueEmu_SelectedETrueNorm"));
  TH2D *ETrueERecNeutron_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueERecNeutron_Selected->Clone("ETrueERecNeutron_SelectedETrueNorm"));
  TH2D *ETrueERecNeutron_timesep_Selected_ETrueNorm =
      static_cast<TH2D *>(ETrueERecNeutron_timesep_Selected->Clone(
          "ETrueERecNeutron_timesep_SelectedETrueNorm"));
  TH2D *ETrueERecProton_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueERecProton_Selected->Clone("ETrueERecProton_SelectedETrueNorm"));
  TH2D *ETrueERecProton_timesep_Selected_ETrueNorm =
      static_cast<TH2D *>(ETrueERecProton_timesep_Selected->Clone(
          "ETrueERecProton_timesep_SelectedETrueNorm"));
  TH2D *ETrueERecPion_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueERecPion_Selected->Clone("ETrueERecPion_SelectedETrueNorm"));
  TH2D *ETrueERecPion_timesep_Selected_ETrueNorm =
      static_cast<TH2D *>(ETrueERecPion_timesep_Selected->Clone(
          "ETrueERecPion_timesep_SelectedETrueNorm"));
  TH2D *ETrueERecOther_Selected_ETrueNorm = static_cast<TH2D *>(
      ETrueERecOther_Selected->Clone("ETrueERecOther_SelectedETrueNorm"));
  TH2D *ETrueERecOther_timesep_Selected_ETrueNorm =
      static_cast<TH2D *>(ETrueERecOther_timesep_Selected->Clone(
          "ETrueERecOther_timesep_SelectedETrueNorm"));

  for (Int_t bin_it = 0; bin_it < ETrue_Selected->GetXaxis()->GetNbins();
       ++bin_it) {
    double NETrueSlice = ETrue_Selected->GetBinContent(bin_it + 1);
    if (!NETrueSlice) {
      continue;
    }
    for (Int_t yit = 0;
         yit < ETrueERec_Selected_ETrueNorm->GetYaxis()->GetNbins(); ++yit) {
      EProxyERec_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          EProxyERec_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueEProxy_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueEProxy_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERec_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERec_Selected->GetBinContent(bin_it + 1, yit + 1) / NETrueSlice);
      ETrueEmu_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueEmu_Selected->GetBinContent(bin_it + 1, yit + 1) / NETrueSlice);
      ETrueERecNeutron_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecNeutron_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERecNeutron_timesep_Selected->SetBinContent(
          bin_it + 1, yit + 1, ETrueERecNeutron_timesep_Selected->GetBinContent(
                                   bin_it + 1, yit + 1) /
                                   NETrueSlice);
      ETrueERecProton_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecProton_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERecProton_timesep_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecProton_timesep_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERecPion_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecPion_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERecPion_timesep_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecPion_timesep_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERecOther_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecOther_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
      ETrueERecOther_timesep_Selected->SetBinContent(
          bin_it + 1, yit + 1,
          ETrueERecOther_timesep_Selected->GetBinContent(bin_it + 1, yit + 1) /
              NETrueSlice);
    }
  }

  EProxyERec_Selected_ETrueNorm->SetDirectory(of);
  ETrueEProxy_Selected_ETrueNorm->SetDirectory(of);
  ETrueERec_Selected_ETrueNorm->SetDirectory(of);
  ETrueEmu_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecNeutron_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecNeutron_timesep_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecProton_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecProton_timesep_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecPion_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecPion_timesep_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecOther_Selected_ETrueNorm->SetDirectory(of);
  ETrueERecOther_timesep_Selected_ETrueNorm->SetDirectory(of);

  of->Write();
}
