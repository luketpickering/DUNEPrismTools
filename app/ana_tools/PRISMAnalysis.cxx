#include "DepositsSummaryTreeReader.hxx"
#include "FluxFitResultsTreeReader.hxx"
#include "OscillationHelper.hxx"
#include "SelectionSummaryTreeReader.hxx"
#include "SliceConfigTreeReader.hxx"
#include "StopConfigTreeReader.hxx"

#include "CovarianceHelper.hxx"

#include "GetUsage.hxx"
#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

#include "TMath.h"
#include "TTree.h"

#include <cmath>
#include <string>
#include <vector>
#include <memory>

std::string NDEventFile;
std::string NDFluxThrowFile;
std::string NDXSecThrowFile;
std::string NDEffFile;
std::string FluxFitFile;
double NDInputPOTPerFile = 5E16;
double NDInputPOT;

std::string NDDataFile;

std::string FDEventFile;
std::string FDFluxThrowFile;
std::string FDXSecThrowFile;
std::string FDEffFile;
double FDInputPOTPerFile = 5E20;
double FDInputPOT;

std::string FDDataFile;

std::string OutputFile = "PRISMAnalysis.root";

std::vector<double> ProjectionBinning;

DepositsSummary::ProjectionVar projVar = DepositsSummary::kERec;

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
    } else if (std::string(argv[opt]) == "-F") {
      FluxFitFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-NI") {
      NDEventFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-NF") {
      NDFluxThrowFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-NX") {
      NDXSecThrowFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-NE") {
      NDEffFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-ND") {
      NDDataFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-FI") {
      FDEventFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-FF") {
      FDFluxThrowFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-FX") {
      FDXSecThrowFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-FE") {
      FDEffFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-FD") {
      FDDataFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      OutputFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-b") {
      ProjectionBinning = BuildBinEdges(argv[++opt]);
    } else if (std::string(argv[opt]) == "-PV") {
      switch (str2T<int>(argv[++opt])) {
      case 1: {
        projVar = DepositsSummary::kETrue;
        break;
      }
      case 2: {
        projVar = DepositsSummary::kEAvail_True;
        break;
      }
      case 3: {
        projVar = DepositsSummary::kERec;
        break;
      }
      case 4: {
        projVar = DepositsSummary::kEFSLep_True;
        break;
      }
      default: {
        std::cout << "[ERROR]: Unknown option for -PV: " << argv[opt]
                  << std::endl;
        SayUsage(argv);
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

enum histlistenum { kFDObs = 0, kFluxcorr, kAcccorr, kLinCombND };
void ScalePRISMPredictions(std::vector<TH1D *> const &inputs, bool PeakScale,
                           TDirectory *dir) {

  TH1D *FDObservation_scale = static_cast<TH1D *>(
      inputs[kFDObs]->Clone((std::string(inputs[kFDObs]->GetName()) +
                             (PeakScale ? "_peak_norm" : "_area_norm"))
                                .c_str()));

  TH1D *FluxCorr_scale = static_cast<TH1D *>(
      inputs[kFluxcorr]->Clone((std::string(inputs[kFluxcorr]->GetName()) +
                                (PeakScale ? "_peak_norm" : "_area_norm"))
                                   .c_str()));
  TH1D *AccCorr_scale = static_cast<TH1D *>(
      inputs[kAcccorr]->Clone((std::string(inputs[kAcccorr]->GetName()) +
                               (PeakScale ? "_peak_norm" : "_area_norm"))
                                  .c_str()));

  TH1D *LinCombND_scale = static_cast<TH1D *>(
      inputs[kLinCombND]->Clone((std::string(inputs[kLinCombND]->GetName()) +
                                 (PeakScale ? "_peak_norm" : "_area_norm"))
                                    .c_str()));

  TH1D *ScaleHelper = static_cast<TH1D *>(
      inputs[kFDObs]->Clone((std::string(inputs[kFDObs]->GetName()) +
                             (PeakScale ? "_peak_norm" : "_area_norm"))
                                .c_str()));
  ScaleHelper->SetDirectory(nullptr);
  ScaleHelper->Add(FluxCorr_scale, -1);
  ScaleHelper->Add(AccCorr_scale, -1);

  double FDToNDScale =
      PeakScale ? (LinCombND_scale->GetMaximum() / ScaleHelper->GetMaximum())
                : (LinCombND_scale->Integral() / ScaleHelper->Integral());

  FDObservation_scale->Scale(FDToNDScale);
  FluxCorr_scale->Scale(FDToNDScale);
  AccCorr_scale->Scale(FDToNDScale);

  TH1D *LinCombND_WithCorr = static_cast<TH1D *>(inputs[kLinCombND]->Clone(
      (std::string(inputs[kLinCombND]->GetName()) +
       (PeakScale ? "_with_corr_peak_norm" : "_with_corr_area_norm"))
          .c_str()));
  LinCombND_WithCorr->Add(FluxCorr_scale);
  LinCombND_WithCorr->Add(AccCorr_scale);

  FDObservation_scale->SetDirectory(dir);
  FluxCorr_scale->SetDirectory(dir);
  AccCorr_scale->SetDirectory(dir);
  LinCombND_scale->SetDirectory(dir);
  LinCombND_WithCorr->SetDirectory(dir);
}

TH1D *GetFDFluxCorrectionHelper() {
  TH1D *BestFit = GetHistogram<TH1D>(FluxFitFile, "BestFit");
  TH1D *InputFlux = GetHistogram<TH1D>(FluxFitFile, "InputFlux");

  TH1D *InputBestFitRatio =
      static_cast<TH1D *>(BestFit->Clone("InputBestFitRatio"));
  InputBestFitRatio->Reset();
  InputBestFitRatio->SetDirectory(nullptr);

  for (Int_t bi_it = 0; bi_it < InputBestFitRatio->GetXaxis()->GetNbins();
       ++bi_it) {
    if (!InputFlux->GetBinContent(bi_it + 1)) {
      continue;
    }
    double frac = 1.0 - (BestFit->GetBinContent(bi_it + 1) /
                         InputFlux->GetBinContent(bi_it + 1));
    InputBestFitRatio->SetBinContent(bi_it + 1, frac);
  }

  delete BestFit;
  delete InputFlux;

  return InputBestFitRatio;
}

TH1D *GetOverlapWeightHelper(Double_t const (&VertexSelectionFV)[3]) {

  StopConfig sc(NDEventFile);
  sc.DetermineNStops();

  std::vector<BoundingBox> StopActiveRegions = sc.GetStopBoundingBoxes(
      true, {VertexSelectionFV[0], VertexSelectionFV[1], VertexSelectionFV[2]});

  Double_t MaxAbsX = -std::numeric_limits<double>::max();
  Double_t MinAbsX = std::numeric_limits<double>::max();

  for (BoundingBox const &bb : StopActiveRegions) {
    MaxAbsX = std::max(MaxAbsX, bb.Max[0]);
    MinAbsX = std::min(MinAbsX, bb.Min[0]);
  }

  size_t NSteps = lrint(MaxAbsX - MinAbsX) / 2;

  std::cout << "[INFO]: Checking coverage [ " << MinAbsX << " -- " << MaxAbsX
            << " : " << ((MaxAbsX - MinAbsX) / double(NSteps))
            << " ], NSteps = " << NSteps << std::endl;

  TH1D *XRangePOTExposureHelper =
      new TH1D("XRangePOTExposureHelper", "", NSteps, MinAbsX, MaxAbsX);
  XRangePOTExposureHelper->SetDirectory(nullptr);

  for (BoundingBox const &bb : StopActiveRegions) {

    for (Int_t bi_it = 0;
         bi_it < XRangePOTExposureHelper->GetXaxis()->GetNbins(); ++bi_it) {

      bool BinCenterIn =
          ((XRangePOTExposureHelper->GetXaxis()->GetBinCenter(bi_it + 1) >
            bb.Min[0]) &&
           (XRangePOTExposureHelper->GetXaxis()->GetBinCenter(bi_it + 1) <
            bb.Max[0]));

      if (BinCenterIn) {
        XRangePOTExposureHelper->AddBinContent(bi_it + 1, sc.POTExposure);
      }
    }
  }

  sc.GetEntry(0);
  for (Int_t bi_it = 0; bi_it < XRangePOTExposureHelper->GetXaxis()->GetNbins();
       ++bi_it) {
    if (!XRangePOTExposureHelper->GetBinContent(bi_it + 1)) {
      continue;
    }

    XRangePOTExposureHelper->SetBinContent(
        bi_it + 1,
        sc.POTExposure / XRangePOTExposureHelper->GetBinContent(bi_it + 1));
  }

  return XRangePOTExposureHelper;
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);
  if (!NDEventFile.size() || !FluxFitFile.size() || !FDEventFile.size()) {
    std::cout
        << "[ERROR]: Was not passed an ND and FD input analysis file and a "
           "flux fit file."
        << std::endl;
  }

  if (!ProjectionBinning.size()) {
    int argc_dum = 3;
    char const *argv_dum[] = {"", "-b", "0_10:0.25"};
    handleOpts(argc_dum, argv_dum);
  }

  size_t MaxThrow_ND, MinThrow_ND;
  double MaxThrow_ND_val = -std::numeric_limits<double>::max(),
         MinThrow_ND_val = std::numeric_limits<double>::max();

  // Build Flux slice helper
  SliceConfig slCfg(FluxFitFile);

  std::vector<std::pair<double, double>> XRange_comp = slCfg.GetXRanges();
  std::vector<double> Coeffs_comp = slCfg.GetCoeffs();

  std::pair<std::vector<double>, std::vector<double>> xrb =
      SliceConfig::BuildXRangeBinsCoeffs(XRange_comp, Coeffs_comp.data(), true);

  std::vector<double> &XRangeBins = xrb.first;
  std::vector<double> &Coeffs = xrb.second;

  if (!XRangeBins.size()) {
    std::cout << "[ERROR]: Got no XRangeBins." << std::endl;
    throw;
  }

  TH1D *CoeffWeightingHelper = new TH1D(
      "CoeffWeightingHelper", "", (XRangeBins.size() - 1), XRangeBins.data());
  CoeffWeightingHelper->SetDirectory(nullptr);
  for (size_t bin_it = 1; bin_it < XRangeBins.size(); ++bin_it) {
    CoeffWeightingHelper->SetBinContent(bin_it, Coeffs[bin_it - 1]);
  }

  // Build FD oscillation weights
  DepositsSummary FDEventReader(FDEventFile);
  SelectionSummary FDSelection(FDEventFile);
  Long64_t NEntries = FDEventReader.GetEntries();

  OscillationHelper FDOscHelper;
  FDOscHelper.Setup(FluxFitFile);

  // Build FD flux correction helper
  TH1D *FDFluxComponentHelper = GetFDFluxCorrectionHelper();

  // ******** Set Up FD Friend Trees
  // Set up Eff tree
  TChain *FDEffFriendTree = nullptr;
  double FDEffWeight = 1;
  if (FDEffFile.size()) {
    FDEffFriendTree = OpenTChainWithFileList("EffWeights", FDEffFile);
    if (FDEffFriendTree) {
      FDEffFriendTree->SetBranchAddress("EffWeight", &FDEffWeight);
    }
  }
  // Set up systematics
  Int_t FDNSystThrows = 0;
  // Set up XSecThrows tree
  TChain *FDXSecThrowsTree = nullptr;
  Int_t FDNXSecThrows = 0;
  double *FDXSecWeights = nullptr;
  if (FDXSecThrowFile.size()) {
    TTree *FDXSecConfigTree =
        OpenTChainWithFileList("ConfigTree", FDXSecThrowFile);
    if (FDXSecConfigTree) {
      FDXSecConfigTree->SetBranchAddress("NThrows", &FDNXSecThrows);
      FDXSecConfigTree->GetEntry(0);
      std::cout << "[INFO]: Read FD XSec tree with " << FDNXSecThrows
                << " throws." << std::endl;
      delete FDXSecConfigTree;

      FDXSecThrowsTree = OpenTChainWithFileList("XSecWeights", FDXSecThrowFile);
      if (!FDXSecThrowsTree) {
        exit(1);
      }

      FDXSecWeights = new double[FDNXSecThrows];
      FDXSecThrowsTree->SetBranchAddress("XSecWeights", FDXSecWeights);
    }
  }

  // Set up FluxThrows
  TChain *FDFluxThrowsTree = nullptr;
  Int_t FDNFluxThrows = std::numeric_limits<Int_t>::max();
  double *FDFluxWeights = nullptr;

  // Set up Thrown distributions
  FDNSystThrows = std::min(FDNXSecThrows, FDNFluxThrows);
  std::vector<TH1D *> ThrownFDObservations;
  std::vector<TH1D *> ThrownFluxCorrections;
  for (Int_t t_it = 0; t_it < FDNSystThrows; ++t_it) {
    ThrownFDObservations.push_back(new TH1D(
        (std::string("FDObservation_") + to_str(projVar) + "_t").c_str(),
        (std::string(";") + to_title(projVar) + ";Count").c_str(),
        (ProjectionBinning.size() - 1), ProjectionBinning.data()));

    ThrownFDObservations.back()->SetDirectory(nullptr);

    ThrownFluxCorrections.push_back(new TH1D(
        (std::string("FDFluxCorrection_") + to_str(projVar) + "_t").c_str(),
        (std::string(";") + to_title(projVar) + ";Count").c_str(),
        (ProjectionBinning.size() - 1), ProjectionBinning.data()));

    ThrownFluxCorrections.back()->SetDirectory(nullptr);
  }

  // Get the acceptance cut used in the ND selection to build the corrector.
  SelectionSummary NDSelection(NDEventFile);

  // Build nominal FD prediction
  TH1D *FDObservation =
      new TH1D((std::string("FDObservation_") + to_str(projVar)).c_str(),
               (std::string(";") + to_title(projVar) + ";Count").c_str(),
               (ProjectionBinning.size() - 1), ProjectionBinning.data());
  FDObservation->SetDirectory(nullptr);

  TH1D *FDFluxCorrection =
      new TH1D((std::string("FDFluxCorrection_") + to_str(projVar)).c_str(),
               (std::string(";") + to_title(projVar) + ";Count").c_str(),
               (ProjectionBinning.size() - 1), ProjectionBinning.data());
  FDFluxCorrection->SetDirectory(nullptr);

  TH1D *FDAcceptanceCorrection = new TH1D(
      (std::string("FDAcceptanceCorrection_") + to_str(projVar)).c_str(),
      (std::string(";") + to_title(projVar) + ";Count").c_str(),
      (ProjectionBinning.size() - 1), ProjectionBinning.data());
  FDAcceptanceCorrection->SetDirectory(nullptr);

  double EHadrVis_AcceptanceCut_GeV =
      NDSelection.EHadrVis_AcceptanceCut_GeV > 1E-5
          ? NDSelection.EHadrVis_AcceptanceCut_GeV
          : std::numeric_limits<double>::max();

  Long64_t LoudEvery = NEntries / 10;
  LoudEvery = LoudEvery ? LoudEvery : 1;
  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    FDEventReader.GetEntry(e_it);
    if (FDEffFriendTree) {
      FDEffFriendTree->GetEntry(e_it);
    }
    if (LoudEvery && !(e_it % LoudEvery)) {
      std::cout << "[INFO]: Processed " << e_it << "/" << NEntries
                << " FD events." << std::endl;
    }

    double FDOscWeight = FDOscHelper.GetWeight(
        FDEventReader.GetProjection(DepositsSummary::kETrue));
    Int_t FluxCorrectionBin = FDFluxComponentHelper->GetXaxis()->FindFixBin(
        FDEventReader.GetProjection(DepositsSummary::kETrue));
    double FluxCorrectionWeight =
        FDFluxComponentHelper->GetBinContent(FluxCorrectionBin);

    FDObservation->Fill(FDEventReader.GetProjection(projVar),
                        FDOscWeight * FDEffWeight);

    FDFluxCorrection->Fill(FDEventReader.GetProjection(projVar),
                           FDOscWeight * FluxCorrectionWeight * FDEffWeight);

    if (FDEventReader.GetProjection(DepositsSummary::kEHadr_vis) >
        EHadrVis_AcceptanceCut_GeV) {
      FDAcceptanceCorrection->Fill(FDEventReader.GetProjection(projVar),
                                   FDOscWeight * FDEffWeight);
    }

    if (FDNSystThrows) {
      if (FDXSecThrowsTree) {
        FDXSecThrowsTree->GetEntry(e_it);
      }
      if (FDFluxThrowsTree) {
        FDFluxThrowsTree->GetEntry(e_it);
      }
      for (Int_t t_it = 0; t_it < FDNSystThrows; ++t_it) {
        double W = 1;
        if (FDXSecThrowsTree) {
          W *= FDXSecWeights[t_it];

          ThrownFluxCorrections[t_it]->Fill(
              FDEventReader.GetProjection(projVar),
              FDOscWeight * FluxCorrectionWeight * FDEffWeight *
                  FDXSecWeights[t_it]);
        }
        if (FDFluxThrowsTree) {
          W *= FDFluxWeights[t_it];
        }

        if (W < 0 || W > 5) {
          std::cout << "[INFO]: FD Throw " << t_it << " event " << e_it
                    << " saw thrown weight " << W << std::endl;
          if (FDXSecThrowsTree) {
            std::cout << "\t\t XSec: " << FDXSecWeights[t_it] << std::endl;
          }
          if (FDFluxThrowsTree) {
            std::cout << "\t\t Flux: " << FDFluxWeights[t_it] << std::endl;
          }
        }

        ThrownFDObservations[t_it]->Fill(FDEventReader.GetProjection(projVar),
                                         FDOscWeight * FDEffWeight * W);
      }
    }
  }

  // ******** Set Up ND Friend Trees
  // Set up Eff tree
  TChain *NDEffFriendTree = nullptr;
  double NDEffWeight = 1;
  if (NDEffFile.size()) {
    NDEffFriendTree = OpenTChainWithFileList("EffWeights", NDEffFile);
    if (NDEffFriendTree) {
      NDEffFriendTree->SetBranchAddress("EffWeight", &NDEffWeight);
    }
  }
  // Set up systematics
  Int_t NSystThrows = 0;
  // Set up XSecThrows tree
  TChain *NDXSecThrowsTree = nullptr;
  Int_t NDNXSecThrows = 0;
  double *NDXSecWeights = nullptr;
  if (NDXSecThrowFile.size()) {
    TTree *NDXSecConfigTree =
        OpenTChainWithFileList("ConfigTree", NDXSecThrowFile);
    if (NDXSecConfigTree) {
      NDXSecConfigTree->SetBranchAddress("NThrows", &NDNXSecThrows);
      NDXSecConfigTree->GetEntry(0);
      std::cout << "[INFO]: Read ND XSec tree with " << NDNXSecThrows
                << " throws." << std::endl;
      delete NDXSecConfigTree;

      NDXSecThrowsTree = OpenTChainWithFileList("XSecWeights", NDXSecThrowFile);
      if (!NDXSecThrowsTree) {
        exit(1);
      }

      NDXSecWeights = new double[NDNXSecThrows];
      NDXSecThrowsTree->SetBranchAddress("XSecWeights", NDXSecWeights);
    }
  }

  // Set up FluxThrows
  TChain *NDFluxThrowsTree = nullptr;
  Int_t NDNFluxThrows = std::numeric_limits<Int_t>::max();
  double *NDFluxWeights = nullptr;

  // Set up Thrown distributions
  NSystThrows = std::min(NDNXSecThrows, std::min(NDNFluxThrows, FDNSystThrows));

  if (NSystThrows) {
    std::cout << "[INFO]: Processing " << NSystThrows
              << " systematic weight throws." << std::endl;
  }

  std::vector<TH1D *> ThrownNDPredictions;
  for (Int_t t_it = 0; t_it < NSystThrows; ++t_it) {
    ThrownNDPredictions.push_back(new TH1D(
        (std::string("NDObservation_") + to_str(projVar) + "_t").c_str(),
        (std::string(";") + to_title(projVar) + ";Count").c_str(),
        (ProjectionBinning.size() - 1), ProjectionBinning.data()));

    ThrownNDPredictions.back()->SetDirectory(nullptr);
  }

  TH1D *LinCombNDObservation =
      new TH1D((std::string("LinCombNDObservation_") + to_str(projVar)).c_str(),
               (std::string(";") + to_title(projVar) + ";Count").c_str(),
               (ProjectionBinning.size() - 1), ProjectionBinning.data());
  LinCombNDObservation->SetDirectory(nullptr);

  DepositsSummary NDEventReader(NDEventFile);

  // ****************** START building ND prediction.
  // Build Overlap weighter
  TH1D *OverlapWeightHelper =
      GetOverlapWeightHelper(NDSelection.VertexSelectionFV);

  LoudEvery = NDEventReader.GetEntries() / 10;
  NEntries = NDEventReader.GetEntries();

  std::vector<TH1D *> Slice;
  std::vector<TH1D *> Slice_Weighted;
  for (size_t slice_it = 0; slice_it < (XRangeBins.size() - 1); ++slice_it) {
    Slice.push_back(new TH1D(
        (std::string("Slice") + to_str(slice_it) + "_" + to_str(projVar))
            .c_str(),
        (std::string(";") + to_title(projVar) + ";Count").c_str(),
        (ProjectionBinning.size() - 1), ProjectionBinning.data()));
    Slice_Weighted.push_back(
        new TH1D((std::string("Slice") + to_str(slice_it) + "_" +
                  to_str(projVar) + "_LCWeighted")
                     .c_str(),
                 (std::string(";") + to_title(projVar) + ";Count").c_str(),
                 (ProjectionBinning.size() - 1), ProjectionBinning.data()));

    Slice.back()->SetDirectory(nullptr);
    Slice_Weighted.back()->SetDirectory(nullptr);
  }

  NDInputPOT = NDSelection.TotalPOT;
  FDInputPOT = FDSelection.TotalPOT;

  Long64_t Fills = 0;
  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    NDEventReader.GetEntry(e_it);

    if (LoudEvery && !(e_it % LoudEvery)) {
      std::cout << "[INFO]: Processed " << e_it << "/" << NEntries
                << " ND events." << std::endl;
    }

    Int_t xb =
        CoeffWeightingHelper->GetXaxis()->FindFixBin(NDEventReader.vtx[0]);

    double CoeffWeight = CoeffWeightingHelper->GetBinContent(xb);

    if (!CoeffWeight) {
      continue;
    }

    if (NDEffFriendTree) {
      NDEffFriendTree->GetEntry(e_it);
    }

    Fills++;

    Int_t xrb =
        OverlapWeightHelper->GetXaxis()->FindFixBin(NDEventReader.vtx[0]);

    double NDOverlapWeight = OverlapWeightHelper->GetBinContent(xrb);
    double NDNominalEventWeight =
        NDEventReader.stop_weight * NDEffWeight * NDOverlapWeight * CoeffWeight;

    Slice[xb - 1]->Fill(NDEventReader.GetProjection(projVar),
                        NDEventReader.stop_weight * NDEffWeight *
                            NDOverlapWeight);
    Slice_Weighted[xb - 1]->Fill(NDEventReader.GetProjection(projVar),
                                 NDNominalEventWeight);

    LinCombNDObservation->Fill(NDEventReader.GetProjection(projVar),
                               NDNominalEventWeight);

    if (NSystThrows) {
      if (NDXSecThrowsTree) {
        NDXSecThrowsTree->GetEntry(e_it);
      }
      if (NDFluxThrowsTree) {
        NDFluxThrowsTree->GetEntry(e_it);
      }
      for (Int_t t_it = 0; t_it < NSystThrows; ++t_it) {
        double W = 1;
        if (NDXSecThrowsTree) {
          W *= NDXSecWeights[t_it];
        }
        if (NDFluxThrowsTree) {
          W *= NDFluxWeights[t_it];
        }

        if (W < 0 || W > 5) {
          std::cout << "[INFO]: ND Throw " << t_it << " event " << e_it
                    << " saw thrown weight " << W << std::endl;
          if (NDXSecThrowsTree) {
            std::cout << "\t\t XSec: " << NDXSecWeights[t_it] << std::endl;
          }
          if (NDFluxThrowsTree) {
            std::cout << "\t\t Flux: " << NDFluxWeights[t_it] << std::endl;
          }
        }

        ThrownNDPredictions[t_it]->Fill(NDEventReader.GetProjection(projVar),
                                        NDNominalEventWeight * W);
      }
    }
  }

  if (NSystThrows) {
    for (Int_t t_it = 0; t_it < NSystThrows; ++t_it) {
      if (ThrownNDPredictions[t_it]->GetMaximum() > MaxThrow_ND_val) {
        MaxThrow_ND = t_it;
        MaxThrow_ND_val = ThrownNDPredictions[t_it]->GetMaximum();
      }
      if (ThrownNDPredictions[t_it]->GetMaximum() < MinThrow_ND_val) {
        MinThrow_ND = t_it;
        MinThrow_ND_val = ThrownNDPredictions[t_it]->GetMaximum();
      }
    }
  }

  std::cout << "[INFO]: Predictions built out of " << Fills << " fills."
            << std::endl;

  TFile *outfile = CheckOpenFile(OutputFile, "RECREATE");

  CoeffWeightingHelper->SetDirectory(outfile);
  FDObservation->SetDirectory(outfile);
  FDFluxCorrection->SetDirectory(outfile);
  FDAcceptanceCorrection->SetDirectory(outfile);
  OverlapWeightHelper->SetDirectory(outfile);
  LinCombNDObservation->SetDirectory(outfile);

  TDirectory *oupD = outfile;

  TDirectory *sliceD = oupD->mkdir("SliceDistributions");
  sliceD->cd();

  for (size_t slice_it = 0; slice_it < (XRangeBins.size() - 1); ++slice_it) {
    Slice[slice_it]->SetDirectory(sliceD);
    Slice_Weighted[slice_it]->SetDirectory(sliceD);
  }

  oupD->cd();

  TDirectory *resD = oupD->mkdir("Results");
  resD->cd();

  if (NSystThrows) {
    CovarianceBuilder NDFDFreedomCovMat(ProjectionBinning.size() - 1);

    std::vector<std::unique_ptr<TH1D>> NDFDDifferences;
    for (Int_t t_it = 0; t_it < NSystThrows; ++t_it) {
      TH1D *NDFDDifference =
          static_cast<TH1D *>(ThrownNDPredictions.front()->Clone());
      NDFDDifference->Reset();
      NDFDDifference->SetDirectory(nullptr);

      // Remove ND corrections built from FD MC from FD 'prediction'
      NDFDDifference->Add(ThrownFDObservations[t_it],
                          ThrownFluxCorrections[t_it], 1, -1);
      NDFDDifference->Add(FDAcceptanceCorrection, -1);

      std::cout << "[INFO]: Throw scale: FDThrowCorr: "
                << NDFDDifference->GetMaximum()
                << ", NDThrow: " << ThrownNDPredictions[t_it]->GetMaximum()
                << std::endl;

      // Get re-scale
      double peak_scale = ThrownNDPredictions[t_it]->GetMaximum() /
                          NDFDDifference->GetMaximum();

      // Re-scale FD to ND and sum
      NDFDDifference->Add(ThrownNDPredictions[t_it], NDFDDifference, 1,
                          -1 * peak_scale);

      NDFDFreedomCovMat.AddThrow_MeanCalc(NDFDDifference);
      NDFDDifferences.emplace_back(NDFDDifference);
    }

    NDFDFreedomCovMat.FinalizeMeanCalc();

    for (Int_t t_it = 0; t_it < NSystThrows; ++t_it) {
      NDFDFreedomCovMat.AddThrow_CovMatCalc(NDFDDifferences[t_it].get());
    }

    NDFDFreedomCovMat.FinalizeCovMatCalc();
    NDFDFreedomCovMat.Write();
  }

  // StopConfig NDSC(NDEventFile);
  // StopConfig FDSC(FDEventFile);
  //
  // double FDFV_m3 = 1;
  // double NDFV_m3 = 1;
  //
  // for (size_t dim_it = 0; dim_it < 3; ++dim_it) {
  //   double NDDim = (NDSC.ActiveMax[dim_it] - NDSC.ActiveMin[dim_it] -
  //                   2.0 * NDSC.VetoGap[dim_it] -
  //                   2.0 * NDSelection.VertexSelectionFV[dim_it]);
  //   double FDDim = (FDSC.ActiveMax[dim_it] - FDSC.ActiveMin[dim_it] -
  //                   2.0 * FDSC.VetoGap[dim_it] -
  //                   2.0 * FDSelection.VertexSelectionFV[dim_it]);
  //   NDFV_m3 *= NDDim;
  //   FDFV_m3 *= FDDim;
  //   std::cout << "[INFO]: ND FV dim[" << dim_it << "] = " << NDDim <<
  //   std::endl; std::cout << "[INFO]: FD FV dim[" << dim_it << "] = " << FDDim
  //   << std::endl;
  // }
  //
  // FluxFitResultsTreeReader op(FluxFitFile);
  //
  // std::cout << "[INFO]: Using ND/FD flux ratio: " <<
  // op.NDOverFDFitScaleFactor
  //           << std::endl;
  // std::cout << "[INFO]: Using ND/FD volume ratio: " << (NDFV_m3 / FDFV_m3)
  //           << std::endl;
  // double FDToNDScale = (NDFV_m3 / FDFV_m3) * op.NDOverFDFitScaleFactor /
  // FDInputPOT;

  ScalePRISMPredictions({FDObservation, FDFluxCorrection,
                         FDAcceptanceCorrection, LinCombNDObservation},
                        true, resD);
  ScalePRISMPredictions({FDObservation, FDFluxCorrection,
                         FDAcceptanceCorrection, LinCombNDObservation},
                        false, resD);

  if (NSystThrows) {
    TDirectory *MaxNDThrow = oupD->mkdir("MaxNDThrow");
    MaxNDThrow->cd();
    TDirectory *MinNDThrow = oupD->mkdir("MinNDThrow");
    MinNDThrow->cd();

    ScalePRISMPredictions(
        {ThrownFDObservations[MaxThrow_ND], ThrownFluxCorrections[MaxThrow_ND],
         FDAcceptanceCorrection, ThrownNDPredictions[MaxThrow_ND]},
        true, MaxNDThrow);
    ScalePRISMPredictions(
        {ThrownFDObservations[MinThrow_ND], ThrownFluxCorrections[MinThrow_ND],
         FDAcceptanceCorrection, ThrownNDPredictions[MinThrow_ND]},
        true, MinNDThrow);
  }

  outfile->Write();
}
