#include "EDepTreeReader.h"

#include "DetectorStop.hxx"
#include "Utils.hxx"

#include "TMath.h"
#include "TTree.h"

#include <string>
#include <vector>

std::string inpfile;
std::string lincombfile;
std::string fluxthrowfile;
std::string xsecthrowfile;
std::string efffile;
std::string outputfile = "PRISMAnalysis.root";
bool SelectOnLepExit = false;
bool UseEDepWeight = false;
bool UseTrueENu = true;
double hadrveto_value = 10E-3;  // GeV

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-i <fulldetprocess.root>             : TChain descriptor for"
               " input tree. \n"
               "\t-F <fluxfitresults.root>             : Result of "
               "dp_FitFluxes to use to build observation.\n"
               "\t-uF <fluxthrowfriend.root>           : Friend tree for -i "
               "option containing flux throws.\n"
               "\t-uX <xsecthrowfriend.root>           : Friend tree for -i "
               "option containing xsec throws.\n"
               "\t-E <efficiencyfriend.root>           : Friend tree for -i "
               "option containing efficiency correction weights.\n"
               "\t-o <outputfile.root>                 : Output file name.\n"
               "\t-v <veto threshold>                  : Hadronic shower veto "
               "threshold in MeV.\n"
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
    } else if (std::string(argv[opt]) == "-uF") {
      fluxthrowfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-uX") {
      xsecthrowfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-E") {
      efffile = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      outputfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-v") {
      hadrveto_value = str2T<double>(argv[++opt]) * 1E-3;
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

  if (!inpfile.size() || !lincombfile.size()) {
    std::cout << "[ERROR]: Was not passed both an input analysis file and a "
                 "flux fit file."
              << std::endl;
  }

  TChain *CoeffTree = new TChain("CoeffTree");

  CoeffTree->Add(lincombfile.c_str());

  double XRange[2];
  double Coeff;

  CoeffTree->SetBranchAddress("XRange", &XRange);
  CoeffTree->SetBranchAddress("Coeff", &Coeff);

  std::vector<double> XRangeBins;
  std::vector<double> Coeffs;

  std::vector<TH1D *> NDObs, NDObs_OAA;

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
      NDObs.push_back(nullptr);
      NDObs_OAA.push_back(nullptr);
    }

    Coeffs.push_back(Coeff);
    XRangeBins.push_back(XRange[1]);
    std::cout << ", " << XRangeBins.back() << std::flush;

    NDObs.push_back(new TH1D(
        (std::string("Obs_") + to_str(XRange[0]) + "_to_" + to_str(XRange[1]))
            .c_str(),
        ";E_{Rec} = E_{#mu} + E_{Hadr,Veto} (GeV);", 100, 0, 10));

    NDObs_OAA.push_back(new TH1D(
        (std::string("OAA_") + to_str(XRange[0]) + "_to_" + to_str(XRange[1]))
            .c_str(),
        ";Off-axis angle (mrad);", 1000, 0, 70));
  }

  std::cout << std::endl;

  TChain *EffFriendTree = nullptr;
  double DumbEffWeight_FV = 1, PosEffWeight_TrueHadrE = 1;

  if (efffile.size()) {
    EffFriendTree = new TChain("EffWeights");

    EffFriendTree->Add(efffile.c_str());

    EffFriendTree->SetBranchAddress("DumbEffWeight_FV", &DumbEffWeight_FV);
    EffFriendTree->SetBranchAddress("PosEffWeight_TrueHadrE",
                                    &PosEffWeight_TrueHadrE);
  }

  TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");

  TH1D *CoeffWeightingHelper = new TH1D(
      "CoeffWeightingHelper", "", XRangeBins.size() - 1, XRangeBins.data());

  TH1D *EventRates =
      static_cast<TH1D *>(CoeffWeightingHelper->Clone("EventRates"));
  TH1D *EventRates_True =
      static_cast<TH1D *>(CoeffWeightingHelper->Clone("EventRates_True"));

  TH1D *EffEventRates = nullptr;

  if (EffFriendTree) {
    EffEventRates =
        static_cast<TH1D *>(CoeffWeightingHelper->Clone("EffEventRates"));
  }

  TH1D *WEffEventRates =
      static_cast<TH1D *>(CoeffWeightingHelper->Clone("WEffEventRates"));

  for (size_t bin_it = 1; bin_it < XRangeBins.size(); ++bin_it) {
    CoeffWeightingHelper->SetBinContent(bin_it, Coeffs[bin_it - 1]);
  }

  TH2D *LinComps = new TH2D(
      "LinComps", ";E_{Rec} = E_{#mu} + E_{Hadr,Veto} (GeV);Offset (cm)", 100,
      0, 10, XRangeBins.size() - 1, XRangeBins.data());

  TH2D *LinComps_eff = static_cast<TH2D *>(LinComps->Clone("EffLinComps"));

  TH2D *LinComps_true = static_cast<TH2D *>(LinComps->Clone("TrueLinComps"));

  EDep edr("EDeps", inpfile);
  size_t loud_every = edr.GetEntries() / 10;
  Long64_t NEntries = edr.GetEntries();

  Long64_t Fills = 0, TFills = 0;
  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    edr.GetEntry(e_it);
    edr.stop_weight = 1;

    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << edr.vtx[0] << ", " << edr.vtx[1] << ", " << edr.vtx[2]
                << "}, Enu: " << edr.nu_4mom[3] << " )" << std::endl;
    }

    Int_t xb = CoeffWeightingHelper->GetXaxis()->FindFixBin(-1 * edr.vtx[0]);

    if ((edr.stop > -1) && (edr.PrimaryLepPDG == 13)) {
      if (!CoeffWeightingHelper->GetBinContent(xb)) {
        continue;
      }
      LinComps_true->Fill(edr.nu_4mom[3], -1 * edr.vtx[0], edr.stop_weight);
      EventRates_True->Fill(-1 * edr.vtx[0], edr.stop_weight);
      TFills++;

      if ((xb > 0) && (xb < CoeffWeightingHelper->GetXaxis()->GetNbins() + 1)) {
        if (NDObs[xb - 1]) {
          NDObs[xb - 1]->Fill(
              UseTrueENu
                  ? edr.nu_4mom[3]
                  : edr.PrimaryLep_4mom[3] + edr.TotalNonlep_Dep_FV +
                        edr.TotalNonlep_Dep_veto,
              CoeffWeightingHelper->GetBinContent(xb) * edr.stop_weight);
        }
        if (NDObs_OAA[xb - 1]) {
          double theta = atan(fabs(edr.vtx[0]) / 57500);
          NDObs_OAA[xb - 1]->Fill(theta * 1.0E3, edr.stop_weight);
        }
      }
    }

    if ((edr.stop < 0) || (SelectOnLepExit && (!edr.LepExit_AboveThresh)) ||
        (!(edr.TotalNonlep_Dep_veto < hadrveto_value)) ||
        (edr.PrimaryLepPDG != 13)) {
      continue;
    }

    if (!CoeffWeightingHelper->GetBinContent(xb)) {
      continue;
    }

    if (EffFriendTree) {
      EffFriendTree->GetEntry(e_it);
    }

    Fills++;

    LinComps->Fill(
        UseTrueENu ? edr.nu_4mom[3] : edr.PrimaryLep_4mom[3] +
                                          edr.TotalNonlep_Dep_FV +
                                          edr.TotalNonlep_Dep_veto,
        -1 * edr.vtx[0], edr.stop_weight);
    LinComps_eff->Fill(
        UseTrueENu ? edr.nu_4mom[3] : edr.PrimaryLep_4mom[3] +
                                          edr.TotalNonlep_Dep_FV +
                                          edr.TotalNonlep_Dep_veto,
        -1 * edr.vtx[0],
        (UseEDepWeight ? DumbEffWeight_FV : PosEffWeight_TrueHadrE) *
            edr.stop_weight);

    EventRates->Fill(-1 * edr.vtx[0]);
    if (EffFriendTree) {
      EffEventRates->Fill(
          -1 * edr.vtx[0],
          (UseEDepWeight ? DumbEffWeight_FV : PosEffWeight_TrueHadrE) *
              edr.stop_weight);
    }

    WEffEventRates->Fill(
        -1 * edr.vtx[0],
        CoeffWeightingHelper->GetBinContent(xb) *
            (UseEDepWeight ? DumbEffWeight_FV : PosEffWeight_TrueHadrE) *
            edr.stop_weight);
  }

  std::cout << "[INFO]: Predictions built out of " << Fills
            << " fills. (TFills = " << TFills << ")." << std::endl;

  TH1D *NDPrediction =
      new TH1D("LinCombPrediction", ";E;Events / GeV", 100, 0, 10);
  TH1D *NDPrediction_True =
      new TH1D("LinCombPrediction_True", ";E;Events / GeV", 100, 0, 10);
  TH2D *LinComps_Weff =
      static_cast<TH2D *>(LinComps_eff->Clone("WEffLinComps"));

  for (Int_t xbin = 1; xbin < NDPrediction->GetXaxis()->GetNbins() + 1;
       ++xbin) {
    double bc = 0, be = 0, bct = 0, bet = 0;

    for (Int_t ybin = 1; ybin < LinComps_eff->GetYaxis()->GetNbins() + 1;
         ++ybin) {
      bc += LinComps_eff->GetBinContent(xbin, ybin) *
            CoeffWeightingHelper->GetBinContent(ybin);
      LinComps_Weff->SetBinContent(
          xbin, ybin, LinComps_eff->GetBinContent(xbin, ybin) *
                          CoeffWeightingHelper->GetBinContent(ybin));
      // std::cout << "[INFO]: Lin sum ERec = "
      //           << NDPrediction->GetXaxis()->GetBinCenter(xbin) << ", OAP = "
      //           << CoeffWeightingHelper->GetXaxis()->GetBinCenter(ybin)
      //           << " cm, W = " << CoeffWeightingHelper->GetBinContent(ybin)
      //           << ", WNEvs = " << LinComps_Weff->GetBinContent(xbin, ybin)
      //           << " (" << bc << ")." << std::endl;
      be += LinComps_eff->GetBinError(xbin, ybin) *
            CoeffWeightingHelper->GetBinContent(ybin);

      bct += LinComps_true->GetBinContent(xbin, ybin) *
             CoeffWeightingHelper->GetBinContent(ybin);
      bet += LinComps_true->GetBinError(xbin, ybin) *
             CoeffWeightingHelper->GetBinContent(ybin);
    }

    // std::cout << "[INFO]: Lin sum ERec = "
    //           << NDPrediction->GetXaxis()->GetBinCenter(xbin) << " BC = " <<
    //           bc
    //           << "." << std::endl;
    NDPrediction->SetBinContent(xbin, bc);
    NDPrediction->SetBinError(xbin, be);

    NDPrediction_True->SetBinContent(xbin, bct);
    NDPrediction_True->SetBinError(xbin, bet);
  }

  NDPrediction->Scale(1, "width");
  NDPrediction_True->Scale(1, "width");

  TH1D *NDPrediction_unorm =
      static_cast<TH1D *>(NDPrediction->Clone("LinCombPrediction_unorm"));
  NDPrediction_unorm->Scale(1.0 / NDPrediction_unorm->Integral("width"));
  TH1D *NDPrediction_True_unorm = static_cast<TH1D *>(
      NDPrediction_True->Clone("LinCombPrediction_True_unorm"));
  NDPrediction_True_unorm->Scale(1.0 /
                                 NDPrediction_True_unorm->Integral("width"));

  TDirectory *oupD = outfile;
  TDirectory *wD = oupD->mkdir("NDObs");
  wD->cd();

  for (TH1D *ob : NDObs) {
    if (!ob) {
      continue;
    }
    ob->SetDirectory(wD);
    ob->Scale(1, "width");
    TH1D *uc = static_cast<TH1D *>(
        ob->Clone((std::string(ob->GetName()) + "_unitnorm").c_str()));
    uc->SetDirectory(wD);
    uc->Scale(1 / uc->Integral("width"));
  }

  oupD->cd();

  TDirectory *wD_2 = oupD->mkdir("NDObs_OAA");
  wD_2->cd();

  for (TH1D *ob : NDObs_OAA) {
    if (!ob) {
      continue;
    }
    ob->SetDirectory(wD_2);
  }

  outfile->Write();
}
