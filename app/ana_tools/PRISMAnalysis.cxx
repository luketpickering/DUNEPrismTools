#include "DepositsSummaryTreeReader.hxx"
#include "StopConfigTreeReader.hxx"
#include "SliceConfigTreeReader.hxx"
#include "OscillationHelper.hxx"
#include "SelectionSummaryTreeReader.hxx"
#include "FluxFitResultsTreeReader.hxx"

#include "StringParserUtility.hxx"
#include "ROOTUtility.hxx"

#include "TMath.h"
#include "TTree.h"

#include <string>
#include <vector>

std::string NDEventFile;
std::string NDFluxThrowFile;
std::string NDXSecThrowFile;
std::string NDEffFile;
std::string FluxFitFile;
double NDInputPOTPerFile = 5E16;
double NDInputPOT;
double NDHadrVisAcceptance = 0; //GeV

std::string NDDataFile;

std::string FDEventFile;
std::string FDFluxThrowFile;
std::string FDXSecThrowFile;
std::string FDEffFile;
double FDInputPOTPerFile = 5E20;
double FDInputPOT;

std::string FDDataFile;

std::string OutputFile = "PRISMAnalysis.root";

std::vector<double> ERecBinning;

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "\t-F <fluxfitresults.root>        : Result of "
         "dp_FitFluxes to use to build observation.\n"
         "\t-NI <NDEvents.root>            : TChain descriptor for"
         " input tree. \n"
         "\t-NF <NDFluxFile.root>          : Friend tree for -NI "
         "option containing flux throws.\n"
         "\t-NX <NDXSecFile.root>          : Friend tree for -NI "
         "option containing xsec throws.\n"
         "\t-NE <NDEffFile.root>           : Friend tree for -NI "
         "option containing efficiency correction weights.\n"
         "\t-ND <NDDataFile.root>          : File containing ND data "
         "distributions.\n"
         "\t-NA <EHadrVis_GeV>             : Shower acceptance cut in GeV.\n"
         "\t-FI <FDEvents.root>            : TChain descriptor for"
         " FD input tree. \n"
         "\t-FF <FDFluxFile.root>          : Friend tree for -FI "
         "option containing flux throws.\n"
         "\t-FX <FDXSecFile.root>          : Friend tree for -FI "
         "option containing xsec throws.\n"
         "\t-FE <FDEffFile.root>           : Friend tree for -FI "
         "option containing efficiency correction weights.\n"
         "\t-FD <FDDataFile.root>          : File containing FD data "
         "distribution.\n"
         "\t-b <low_up:width[,low_up:width...]>  : ENuBinning descriptor.\n"
         "\t-o <OutputFile.root>                 : Output file name.\n"
      << std::endl;
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
    }else if (std::string(argv[opt]) == "-NA") {
      NDHadrVisAcceptance = str2T<double>(argv[++opt]);
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
      ERecBinning = BuildBinEdges(argv[++opt]);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
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

  std::vector<BoundingBox> StopActiveRegions = sc.GetStopBoundingBoxes(true,
    {VertexSelectionFV[0], VertexSelectionFV[1], VertexSelectionFV[2]});

  Double_t MaxAbsX = -std::numeric_limits<double>::max();
  Double_t MinAbsX = std::numeric_limits<double>::max();

  for(BoundingBox const & bb : StopActiveRegions){
    MaxAbsX = std::max(MaxAbsX, bb.Max[0]);
    MinAbsX = std::min(MinAbsX, bb.Min[0]);
  }

  TH1D *XRangePOTExposureHelper =
      new TH1D("XRangePOTExposureHelper", "", 2000, -3800, 200);
  XRangePOTExposureHelper->SetDirectory(nullptr);

  for(BoundingBox const & bb : StopActiveRegions){

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

  if (!ERecBinning.size()) {
    int argc_dum = 3;
    char const *argv_dum[] = {"", "-b", "0_10:0.25"};
    handleOpts(argc_dum, argv_dum);
  }

  // Build Flux slice helper
  SliceConfig slCfg(FluxFitFile);

  std::vector< std::pair<double, double> > XRange_comp = slCfg.GetXRanges();
  std::vector<double> Coeffs_comp = slCfg.GetCoeffs();

  std::pair< std::vector<double>, std::vector<double> >
    xrb = SliceConfig::BuildXRangeBinsCoeffs(XRange_comp, Coeffs_comp.data(),
      true);

  std::vector<double> &XRangeBins = xrb.first;
  std::vector<double> &Coeffs = xrb.second;

  if(!XRangeBins.size()){
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

  std::vector<double> FDOscWeights;
  OscillationHelper FDOscHelper;
  FDOscHelper.Setup(FluxFitFile);
  FDOscWeights.resize(NEntries);

  // Build FD flux correction helper
  std::vector<double> FluxCorrectionWeights;
  FluxCorrectionWeights.resize(NEntries);
  TH1D *FDFluxComponentHelper = GetFDFluxCorrectionHelper();

  // ******** Set Up FD Friend Trees
  // Set up Eff tree
  TChain *FDEffFriendTree = nullptr;
  double FDEffWeight = 1;
  if (FDXSecThrowFile.size()) {
    FDEffFriendTree = OpenTChainWithFileList("EffWeights", FDXSecThrowFile);
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
      FDXSecThrowsTree->Branch("XSecWeights", FDXSecWeights);
    }
  }

  // Set up FluxThrows
  TChain *FDFluxThrowsTree = nullptr;
  Int_t FDNFluxThrows = std::numeric_limits<Int_t>::max();
  double *FDFluxWeights = nullptr;

  // Set up Thrown distributions
  FDNSystThrows = std::min(FDNXSecThrows, FDNFluxThrows);
  std::vector<TH1D *> ThrownFDPredictions_ERec;
  std::vector<TH1D *> ThrownFDFluxCorrection_ERec;
  for (Int_t t_it = 0; t_it < FDNSystThrows; ++t_it) {
    ThrownFDPredictions_ERec.push_back(
        new TH1D("FDPrediction_ERec_t", ";E_{Dep} (GeV);Count",
                 (ERecBinning.size() - 1), ERecBinning.data()));
    ThrownFDFluxCorrection_ERec.push_back(
        new TH1D("FDFluxCorrection_ERec", ";E_{Dep} (GeV);Count",
                 (ERecBinning.size() - 1), ERecBinning.data()));

    ThrownFDPredictions_ERec.back()->SetDirectory(nullptr);
    ThrownFDFluxCorrection_ERec.back()->SetDirectory(nullptr);
  }

  // Build nominal FD prediction
  TH1D *FDPrediction_ERec =
      new TH1D("FDPrediction_ERec", ";E_{Dep} (GeV);Count",
               (ERecBinning.size() - 1), ERecBinning.data());
  TH1D *FDPrediction_True =
      new TH1D("FDPrediction_True", ";E_{#nu} (GeV);Count",
               (ERecBinning.size() - 1), ERecBinning.data());
  TH1D *FDPrediction_EProxy =
      new TH1D("FDPrediction_EProxy", ";E_{#nu,Proxy} (GeV);Count",
               (ERecBinning.size() - 1), ERecBinning.data());
  FDPrediction_ERec->SetDirectory(nullptr);
  FDPrediction_True->SetDirectory(nullptr);
  FDPrediction_EProxy->SetDirectory(nullptr);

  TH1D *FDFluxCorrection_ERec =
      new TH1D("FDFluxCorrection_ERec", ";E_{Dep} (GeV);Count",
               (ERecBinning.size() - 1), ERecBinning.data());
  TH1D *FDFluxCorrection_True =
      new TH1D("FDFluxCorrection_True", ";E_{#nu} (GeV);Count",
               (ERecBinning.size() - 1), ERecBinning.data());
  TH1D *FDFluxCorrection_EProxy =
      new TH1D("FDFluxCorrection_EProxy", ";E_{#nu,Proxy} (GeV);Count",
               (ERecBinning.size() - 1), ERecBinning.data());
  FDFluxCorrection_ERec->SetDirectory(nullptr);
  FDFluxCorrection_True->SetDirectory(nullptr);
  FDFluxCorrection_EProxy->SetDirectory(nullptr);

  TH1D *FDAcceptanceCorrection_ERec =
      new TH1D("FDAcceptanceCorrection_ERec", ";E_{Dep} (GeV);Count",
               (ERecBinning.size() - 1), ERecBinning.data());
  TH1D *FDAcceptanceCorrection_True =
      new TH1D("FDAcceptanceCorrection_True", ";E_{#nu} (GeV);Count",
               (ERecBinning.size() - 1), ERecBinning.data());
  TH1D *FDAcceptanceCorrection_EProxy =
      new TH1D("FDAcceptanceCorrection_EProxy", ";E_{#nu,Proxy} (GeV);Count",
               (ERecBinning.size() - 1), ERecBinning.data());
  FDAcceptanceCorrection_ERec->SetDirectory(nullptr);
  FDAcceptanceCorrection_True->SetDirectory(nullptr);
  FDAcceptanceCorrection_EProxy->SetDirectory(nullptr);

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

    FDOscWeights[e_it] = FDOscHelper.GetWeight(
      FDEventReader.GetProjection(DepositsSummary::kETrue));
    Int_t FluxCorrectionBin =
        FDFluxComponentHelper->GetXaxis()->FindFixBin(
          FDEventReader.GetProjection(DepositsSummary::kETrue));
    FluxCorrectionWeights[e_it] =
        FDFluxComponentHelper->GetBinContent(FluxCorrectionBin);

    FDPrediction_ERec->Fill(FDEventReader.GetProjection(DepositsSummary::kERec),
                            FDOscWeights[e_it] * FDEffWeight);
    FDPrediction_True->Fill(
      FDEventReader.GetProjection(DepositsSummary::kETrue),
                            FDOscWeights[e_it] * FDEffWeight);
    FDPrediction_EProxy->Fill(
      FDEventReader.GetProjection(DepositsSummary::kEAvail_True),
                              FDOscWeights[e_it] * FDEffWeight);

    FDFluxCorrection_ERec->Fill(
        FDEventReader.GetProjection(DepositsSummary::kERec),
        FDOscWeights[e_it] * FluxCorrectionWeights[e_it] * FDEffWeight);
    FDFluxCorrection_True->Fill(
        FDEventReader.GetProjection(DepositsSummary::kETrue),
        FDOscWeights[e_it] * FluxCorrectionWeights[e_it] * FDEffWeight);
    FDFluxCorrection_EProxy->Fill(
        FDEventReader.GetProjection(DepositsSummary::kEAvail_True),
        FDOscWeights[e_it] * FluxCorrectionWeights[e_it] * FDEffWeight);

    if(FDEventReader.GetProjection(DepositsSummary::kEHadr_vis) >
      NDHadrVisAcceptance){
      FDAcceptanceCorrection_ERec->Fill(
          FDEventReader.GetProjection(DepositsSummary::kERec),
          FDOscWeights[e_it] * FDEffWeight);
      FDAcceptanceCorrection_True->Fill(
          FDEventReader.GetProjection(DepositsSummary::kETrue),
          FDOscWeights[e_it] * FDEffWeight);
      FDAcceptanceCorrection_EProxy->Fill(
          FDEventReader.GetProjection(DepositsSummary::kEAvail_True),
          FDOscWeights[e_it] * FDEffWeight);
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
        }
        if (FDFluxThrowsTree) {
          W *= FDFluxWeights[t_it];
        }

        ThrownFDPredictions_ERec[t_it]->Fill(
            FDEventReader.GetProjection(DepositsSummary::kERec),
            FDOscWeights[e_it] * FDEffWeight * W);
        ThrownFDFluxCorrection_ERec[t_it]->Fill(
            FDEventReader.GetProjection(DepositsSummary::kERec),
            FDOscWeights[e_it] * FluxCorrectionWeights[e_it] * FDEffWeight * W);
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
  Int_t NDNSystThrows = 0;
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
      NDXSecThrowsTree->Branch("XSecWeights", NDXSecWeights);
    }
  }

  // Set up FluxThrows
  TChain *NDFluxThrowsTree = nullptr;
  Int_t NDNFluxThrows = std::numeric_limits<Int_t>::max();
  double *NDFluxWeights = nullptr;

  // Set up Thrown distributions
  NDNSystThrows = std::min(NDNXSecThrows, NDNFluxThrows);
  if (NDNSystThrows != FDNSystThrows) {
    std::cout << "[ERROR]: Inputs contain " << NDNSystThrows
              << " ND syst throws and " << FDNSystThrows
              << " FD syst throws. These must be the same." << std::endl;
    exit(1);
  }

  std::vector<TH1D *> ThrownNDPredictions_ERec;
  for (Int_t t_it = 0; t_it < NDNSystThrows; ++t_it) {
    ThrownNDPredictions_ERec.push_back(
        new TH1D("NDPrediction_ERec_t", ";E_{Dep} (GeV);Count",
                 (ERecBinning.size() - 1), ERecBinning.data()));

    ThrownNDPredictions_ERec.back()->SetDirectory(nullptr);
  }

  TH2D *OffSetVsERec =
      new TH2D("OffSetVsERec",
               ";E_{Rec} = E_{#mu} + E_{Hadr,Veto} (GeV);Offset (cm);Count",
               (ERecBinning.size() - 1), ERecBinning.data(),
               XRangeBins.size() - 1, XRangeBins.data());

  TH2D *OffSetVsEProxy =
      new TH2D("OffSetVsEProxy", ";E_{#nu,Proxy} (GeV);Offset (cm);Count",
               (ERecBinning.size() - 1), ERecBinning.data(),
               XRangeBins.size() - 1, XRangeBins.data());
  OffSetVsERec->SetDirectory(nullptr);
  OffSetVsEProxy->SetDirectory(nullptr);

  TH1D *LinCombOneShotERec =
      new TH1D("LinCombOneShotERec",
               ";E_{Rec} = E_{#mu} + E_{Hadr,Veto} (GeV);Events / GeV",
               (ERecBinning.size() - 1), ERecBinning.data());
  TH1D *LinCombOneShotEProxy =
      new TH1D("LinCombOneShotEProxy", ";E_{#nu,Proxy} (GeV);Events / GeV",
               (ERecBinning.size() - 1), ERecBinning.data());
  LinCombOneShotERec->SetDirectory(nullptr);
  LinCombOneShotEProxy->SetDirectory(nullptr);
  TH1D *LinCombOneShotETrue =
      new TH1D("LinCombOneShotETrue", ";E_{#nu} (GeV);Events / GeV",
               (ERecBinning.size() - 1), ERecBinning.data());
  LinCombOneShotETrue->SetDirectory(nullptr);

  DepositsSummary NDEventReader(NDEventFile);
  SelectionSummary NDSelection(NDEventFile);

  // ****************** START building ND prediction.
  // Build Overlap weighter
  TH1D *OverlapWeightHelper =
    GetOverlapWeightHelper(NDSelection.VertexSelectionFV);

  LoudEvery = NDEventReader.GetEntries() / 10;
  NEntries = NDEventReader.GetEntries();
  std::vector<double> NDOverlapWeights;
  NDOverlapWeights.resize(NEntries);

  std::vector<TH1D *> SliceEProxy;
  std::vector<TH1D *> SliceEProxy_Weighted;
  std::vector<TH1D *> SliceERec;
  std::vector<TH1D *> SliceERec_Weighted;
  for (size_t slice_it = 0; slice_it < (XRangeBins.size() - 1); ++slice_it) {
    SliceEProxy.push_back(
        new TH1D((std::string("EProxy_Slice") + to_str(slice_it)).c_str(),
                 ";E_{#nu,Proxy} (GeV);Count", (ERecBinning.size() - 1),
                 ERecBinning.data()));
    SliceEProxy_Weighted.push_back(new TH1D(
        (std::string("EProxy_WeightedSlice") + to_str(slice_it)).c_str(),
        ";E_{Dep} (GeV);Count", (ERecBinning.size() - 1), ERecBinning.data()));

    SliceEProxy.back()->SetDirectory(nullptr);
    SliceEProxy_Weighted.back()->SetDirectory(nullptr);

    SliceERec.push_back(
        new TH1D((std::string("ERec_Slice") + to_str(slice_it)).c_str(),
                 ";E_{#nu,Proxy} (GeV);Count", (ERecBinning.size() - 1),
                 ERecBinning.data()));
    SliceERec_Weighted.push_back(new TH1D(
        (std::string("ERec_WeightedSlice") + to_str(slice_it)).c_str(),
        ";E_{Dep} (GeV);Count", (ERecBinning.size() - 1), ERecBinning.data()));

    SliceERec.back()->SetDirectory(nullptr);
    SliceERec_Weighted.back()->SetDirectory(nullptr);
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

    NDOverlapWeights[e_it] = OverlapWeightHelper->GetBinContent(xrb);

    OffSetVsERec->Fill(
        NDEventReader.GetProjection(DepositsSummary::kERec),
        NDEventReader.vtx[0],
        NDEventReader.stop_weight * NDEffWeight * NDOverlapWeights[e_it]);
    OffSetVsEProxy->Fill(
        NDEventReader.GetProjection(DepositsSummary::kEAvail_True),
          NDEventReader.vtx[0],
            NDEventReader.stop_weight * NDEffWeight * NDOverlapWeights[e_it]);

    SliceEProxy[xb - 1]->Fill(
        NDEventReader.GetProjection(DepositsSummary::kEAvail_True),
        NDEventReader.stop_weight * NDEffWeight * NDOverlapWeights[e_it]);
    SliceEProxy_Weighted[xb - 1]->Fill(
        NDEventReader.GetProjection(DepositsSummary::kEAvail_True),
          NDEventReader.stop_weight * NDEffWeight * NDOverlapWeights[e_it] * \
            CoeffWeight);

    SliceERec[xb - 1]->Fill(
      NDEventReader.GetProjection(DepositsSummary::kERec),
        NDEventReader.stop_weight * NDEffWeight * NDOverlapWeights[e_it]);
    SliceERec_Weighted[xb - 1]->Fill(
      NDEventReader.GetProjection(DepositsSummary::kERec),
         NDEventReader.stop_weight * NDEffWeight *
             NDOverlapWeights[e_it] * CoeffWeight);

    LinCombOneShotERec->Fill(
      NDEventReader.GetProjection(DepositsSummary::kERec),
       NDEventReader.stop_weight * NDEffWeight * NDOverlapWeights[e_it] * \
       CoeffWeight);
    LinCombOneShotEProxy->Fill(
      NDEventReader.GetProjection(DepositsSummary::kEAvail_True),
        NDEventReader.stop_weight * NDEffWeight *
          NDOverlapWeights[e_it] * CoeffWeight);
    LinCombOneShotETrue->Fill(
      NDEventReader.GetProjection(DepositsSummary::kETrue),
        NDEventReader.stop_weight * NDEffWeight *
            NDOverlapWeights[e_it] * CoeffWeight);

    if (NDNSystThrows) {
      if (NDXSecThrowsTree) {
        NDXSecThrowsTree->GetEntry(e_it);
      }
      if (NDFluxThrowsTree) {
        NDFluxThrowsTree->GetEntry(e_it);
      }
      for (Int_t t_it = 0; t_it < NDNSystThrows; ++t_it) {
        double W = 1;
        if (NDXSecThrowsTree) {
          W *= NDXSecWeights[t_it];
        }
        if (NDFluxThrowsTree) {
          W *= NDFluxWeights[t_it];
        }

        ThrownNDPredictions_ERec[t_it]->Fill(
            NDEventReader.GetProjection(DepositsSummary::kERec),
            NDEventReader.stop_weight * NDEffWeight * NDOverlapWeights[e_it] *
                CoeffWeight * W);
      }
    }
  }

  std::cout << "[INFO]: Predictions built out of " << Fills << " fills."
            << std::endl;

  // if (NDNSystThrows) {
  //   std::vector<TH1D *> NDFDDifferences;
  //   for (Int_t t_it = 0; t_it < NDNSystThrows; ++t_it) {
  //   }
  // }

  TH1D *LinCombERec = new TH1D(
      "LinCombERec", ";E_{Rec} = E_{#mu} + E_{Hadr,Veto} (GeV);Events / GeV",
      (ERecBinning.size() - 1), ERecBinning.data());
  TH1D *LinCombEProxy =
      new TH1D("LinCombEProxy", ";E_{#nu,Proxy} (GeV);Events / GeV",
               (ERecBinning.size() - 1), ERecBinning.data());

  LinCombERec->SetDirectory(nullptr);
  LinCombEProxy->SetDirectory(nullptr);

  for (Int_t xbin = 1; xbin < OffSetVsERec->GetXaxis()->GetNbins() + 1;
       ++xbin) {
    double bc_erec = 0, be_erec = 0, bc_eprox = 0, be_eprox = 0;

    for (Int_t ybin = 1; ybin < OffSetVsERec->GetYaxis()->GetNbins() + 1;
         ++ybin) {
      bc_erec += OffSetVsERec->GetBinContent(xbin, ybin) *
                 CoeffWeightingHelper->GetBinContent(ybin);

      be_erec += OffSetVsERec->GetBinError(xbin, ybin) *
                 CoeffWeightingHelper->GetBinContent(ybin);

      bc_eprox += OffSetVsEProxy->GetBinContent(xbin, ybin) *
                  CoeffWeightingHelper->GetBinContent(ybin);
      be_eprox += OffSetVsEProxy->GetBinError(xbin, ybin) *
                  CoeffWeightingHelper->GetBinContent(ybin);
    }

    LinCombERec->SetBinContent(xbin, bc_erec);
    LinCombERec->SetBinError(xbin, be_erec);

    LinCombEProxy->SetBinContent(xbin, bc_eprox);
    LinCombEProxy->SetBinError(xbin, be_eprox);
  }

  TFile *outfile = CheckOpenFile(OutputFile, "RECREATE");

  CoeffWeightingHelper->SetDirectory(outfile);
  FDPrediction_ERec->SetDirectory(outfile);
  FDPrediction_True->SetDirectory(outfile);
  FDPrediction_EProxy->SetDirectory(outfile);
  FDFluxCorrection_ERec->SetDirectory(outfile);
  FDFluxCorrection_True->SetDirectory(outfile);
  FDFluxCorrection_EProxy->SetDirectory(outfile);
  FDAcceptanceCorrection_ERec->SetDirectory(outfile);
  FDAcceptanceCorrection_True->SetDirectory(outfile);
  FDAcceptanceCorrection_EProxy->SetDirectory(outfile);
  OverlapWeightHelper->SetDirectory(outfile);
  OffSetVsERec->SetDirectory(outfile);
  OffSetVsEProxy->SetDirectory(outfile);
  LinCombERec->SetDirectory(outfile);
  LinCombEProxy->SetDirectory(outfile);
  LinCombOneShotERec->SetDirectory(outfile);
  LinCombOneShotEProxy->SetDirectory(outfile);
  LinCombOneShotETrue->SetDirectory(outfile);

  TDirectory *oupD = outfile;

  TDirectory *sliceD = oupD->mkdir("SliceDistributions");
  sliceD->cd();

  for (size_t slice_it = 0; slice_it < (XRangeBins.size() - 1); ++slice_it) {
    SliceEProxy[slice_it]->SetDirectory(sliceD);
    SliceEProxy_Weighted[slice_it]->SetDirectory(sliceD);
    SliceERec[slice_it]->SetDirectory(sliceD);
    SliceERec_Weighted[slice_it]->SetDirectory(sliceD);
  }

  oupD->cd();

  TDirectory *resD = oupD->mkdir("Results");
  resD->cd();

  StopConfig NDSC(NDEventFile);
  StopConfig FDSC(FDEventFile);

  double FDFV_m3 = 1;
  double NDFV_m3 = 1;

  for(size_t dim_it = 0; dim_it < 3; ++dim_it){
    double NDDim = (NDSC.ActiveMax[dim_it] - NDSC.ActiveMin[dim_it]
      - 2.0*NDSC.VetoGap[dim_it]
      - 2.0*NDSelection.VertexSelectionFV[dim_it]);
    double FDDim = (FDSC.ActiveMax[dim_it] - FDSC.ActiveMin[dim_it]
      - 2.0*FDSC.VetoGap[dim_it]
      - 2.0*FDSelection.VertexSelectionFV[dim_it]);
    NDFV_m3 *= NDDim;
    FDFV_m3 *= FDDim;
    std::cout << "[INFO]: ND FV dim[" << dim_it << "] = " << NDDim << std::endl;
    std::cout << "[INFO]: FD FV dim[" << dim_it << "] = " << FDDim << std::endl;
  }

  FluxFitResultsTreeReader op(FluxFitFile);

  std::cout << "[INFO]: Using ND/FD flux ratio: " << op.NDOverFDFitScaleFactor
            << std::endl;
  std::cout << "[INFO]: Using ND/FD volume ratio: " << (NDFV_m3 / FDFV_m3)
            << std::endl;

  TH1D *FDPrediction_ERec_POTScaled = static_cast<TH1D *>(
      FDPrediction_ERec->Clone("FDPrediction_ERec_POTScaled"));
  FDPrediction_ERec_POTScaled->Scale(
      (NDFV_m3 / FDFV_m3) * op.NDOverFDFitScaleFactor / FDInputPOT, "width");
  FDPrediction_ERec_POTScaled->GetYaxis()->SetTitle("Count per 1 GeV per POT");

  TH1D *FDFluxCorrection_ERec_clone =
      static_cast<TH1D *>(FDFluxCorrection_ERec->Clone());
  FDFluxCorrection_ERec_clone->Scale(
      (NDFV_m3 / FDFV_m3) * op.NDOverFDFitScaleFactor / FDInputPOT, "width");

  TH1D *FDAcceptanceCorrection_ERec_clone =
      static_cast<TH1D *>(FDAcceptanceCorrection_ERec->Clone());
  FDAcceptanceCorrection_ERec_clone->Scale(
      (NDFV_m3 / FDFV_m3) * op.NDOverFDFitScaleFactor / FDInputPOT, "width");

  TH1D *LinCombOneShotERec_POTScaled =
      static_cast<TH1D *>(LinCombOneShotERec->Clone("LinCombOneShotERec_POTScaled"));
  LinCombOneShotERec_POTScaled->Scale(1.0 / NDInputPOT, "width");
  LinCombOneShotERec_POTScaled->GetYaxis()->SetTitle("Count per 1 GeV per POT");

  TH1D *FDPrediction_ERec_POTScaled_HEURSCALE = static_cast<TH1D *>(
      FDPrediction_ERec->Clone("FDPrediction_ERec_POTScaled_HEURSCALE"));

  TH1D *FDFluxCorrection_ERec_HEURSCALE = static_cast<TH1D *>(
      FDFluxCorrection_ERec->Clone("FDFluxCorrection_ERec_HEURSCALE"));
  TH1D *FDAcceptanceCorrection_ERec_HEURSCALE = static_cast<TH1D *>(
      FDAcceptanceCorrection_ERec->Clone("FDAcceptanceCorrection_ERec_HEURSCALE"));
  TH1D *LinCombOneShotERec_POTScaled_HEURSCALE = static_cast<TH1D *>(
      LinCombOneShotERec->Clone("LinCombOneShotERec_POTScaled_HEURSCALE"));

  double FDToNDScale = LinCombOneShotERec_POTScaled_HEURSCALE->GetMaximum() /
                       FDPrediction_ERec_POTScaled_HEURSCALE->GetMaximum();

  FDPrediction_ERec_POTScaled_HEURSCALE->Scale(FDToNDScale);
  FDFluxCorrection_ERec_HEURSCALE->Scale(FDToNDScale);
  FDAcceptanceCorrection_ERec_HEURSCALE->Scale(FDToNDScale);

  TH1D *LinCombOneShotERec_POTScaled_WithCorr_HEURSCALE =
      static_cast<TH1D *>(LinCombOneShotERec_POTScaled_HEURSCALE->Clone(
          "LinCombOneShotERec_POTScaled_WithCorr_HEURSCALE"));

  LinCombOneShotERec_POTScaled_WithCorr_HEURSCALE->Add(
      FDFluxCorrection_ERec_HEURSCALE);
  LinCombOneShotERec_POTScaled_WithCorr_HEURSCALE->Add(
      FDAcceptanceCorrection_ERec_HEURSCALE);

  FDPrediction_ERec_POTScaled->SetDirectory(resD);
  FDFluxCorrection_ERec_clone->SetDirectory(resD);
  FDAcceptanceCorrection_ERec_clone->SetDirectory(resD);
  LinCombOneShotERec_POTScaled->SetDirectory(resD);






  TH1D *FDPrediction_EProxy_POTScaled = static_cast<TH1D *>(
      FDPrediction_EProxy->Clone("FDPrediction_EProxy_POTScaled"));
  FDPrediction_EProxy_POTScaled->Scale(
      (NDFV_m3 / FDFV_m3) * op.NDOverFDFitScaleFactor / FDInputPOT, "width");
  FDPrediction_EProxy_POTScaled->GetYaxis()->SetTitle("Count per 1 GeV per POT");

  TH1D *FDFluxCorrection_EProxy_clone =
      static_cast<TH1D *>(FDFluxCorrection_EProxy->Clone());
  FDFluxCorrection_EProxy_clone->Scale(
      (NDFV_m3 / FDFV_m3) * op.NDOverFDFitScaleFactor / FDInputPOT, "width");

  TH1D *FDAcceptanceCorrection_EProxy_clone =
      static_cast<TH1D *>(FDAcceptanceCorrection_EProxy->Clone());
  FDAcceptanceCorrection_EProxy_clone->Scale(
      (NDFV_m3 / FDFV_m3) * op.NDOverFDFitScaleFactor / FDInputPOT, "width");

  TH1D *LinCombOneShotEProxy_POTScaled =
      static_cast<TH1D *>(LinCombOneShotEProxy->Clone("LinCombOneShotEProxy_POTScaled"));
  LinCombOneShotEProxy_POTScaled->Scale(1.0 / NDInputPOT, "width");
  LinCombOneShotEProxy_POTScaled->GetYaxis()->SetTitle("Count per 1 GeV per POT");

  TH1D *FDPrediction_EProxy_POTScaled_HEURSCALE = static_cast<TH1D *>(
      FDPrediction_EProxy->Clone("FDPrediction_EProxy_POTScaled_HEURSCALE"));

  TH1D *FDFluxCorrection_EProxy_HEURSCALE = static_cast<TH1D *>(
      FDFluxCorrection_EProxy->Clone("FDFluxCorrection_EProxy_HEURSCALE"));
  TH1D *FDAcceptanceCorrection_EProxy_HEURSCALE = static_cast<TH1D *>(
      FDAcceptanceCorrection_EProxy->Clone("FDAcceptanceCorrection_EProxy_HEURSCALE"));
  TH1D *LinCombOneShotEProxy_POTScaled_HEURSCALE = static_cast<TH1D *>(
      LinCombOneShotEProxy->Clone("LinCombOneShotEProxy_POTScaled_HEURSCALE"));

  FDToNDScale = LinCombOneShotEProxy_POTScaled_HEURSCALE->GetMaximum() /
                       FDPrediction_EProxy_POTScaled_HEURSCALE->GetMaximum();

  FDPrediction_EProxy_POTScaled_HEURSCALE->Scale(FDToNDScale);
  FDFluxCorrection_EProxy_HEURSCALE->Scale(FDToNDScale);
  FDAcceptanceCorrection_EProxy_HEURSCALE->Scale(FDToNDScale);

  TH1D *LinCombOneShotEProxy_POTScaled_WithCorr_HEURSCALE =
      static_cast<TH1D *>(LinCombOneShotEProxy_POTScaled_HEURSCALE->Clone(
          "LinCombOneShotEProxy_POTScaled_WithCorr_HEURSCALE"));

  LinCombOneShotEProxy_POTScaled_WithCorr_HEURSCALE->Add(
      FDFluxCorrection_EProxy_HEURSCALE);
  LinCombOneShotEProxy_POTScaled_WithCorr_HEURSCALE->Add(
      FDAcceptanceCorrection_EProxy_HEURSCALE);

  FDPrediction_EProxy_POTScaled->SetDirectory(resD);
  FDFluxCorrection_EProxy_clone->SetDirectory(resD);
  FDAcceptanceCorrection_EProxy_clone->SetDirectory(resD);
  LinCombOneShotEProxy_POTScaled->SetDirectory(resD);






  TH1D *FDPrediction_True_POTScaled = static_cast<TH1D *>(
      FDPrediction_True->Clone("FDPrediction_True_POTScaled"));
  FDPrediction_True_POTScaled->Scale(
      (NDFV_m3 / FDFV_m3) * op.NDOverFDFitScaleFactor / FDInputPOT, "width");
  FDPrediction_True_POTScaled->GetYaxis()->SetTitle("Count per 1 GeV per POT");

  TH1D *FDFluxCorrection_True_clone =
      static_cast<TH1D *>(FDFluxCorrection_True->Clone());
  FDFluxCorrection_True_clone->Scale(
      (NDFV_m3 / FDFV_m3) * op.NDOverFDFitScaleFactor / FDInputPOT, "width");

  TH1D *FDAcceptanceCorrection_True_clone =
      static_cast<TH1D *>(FDAcceptanceCorrection_True->Clone());
  FDAcceptanceCorrection_True_clone->Scale(
      (NDFV_m3 / FDFV_m3) * op.NDOverFDFitScaleFactor / FDInputPOT, "width");

  TH1D *LinCombOneShotETrue_POTScaled =
      static_cast<TH1D *>(LinCombOneShotETrue->Clone("LinCombOneShotETrue_POTScaled"));
  LinCombOneShotETrue_POTScaled->Scale(1.0 / NDInputPOT, "width");
  LinCombOneShotETrue_POTScaled->GetYaxis()->SetTitle("Count per 1 GeV per POT");

  TH1D *FDPrediction_True_POTScaled_HEURSCALE = static_cast<TH1D *>(
      FDPrediction_True->Clone("FDPrediction_True_POTScaled_HEURSCALE"));

  TH1D *FDFluxCorrection_True_HEURSCALE = static_cast<TH1D *>(
      FDFluxCorrection_True->Clone("FDFluxCorrection_True_HEURSCALE"));
  TH1D *FDAcceptanceCorrection_True_HEURSCALE = static_cast<TH1D *>(
      FDAcceptanceCorrection_True->Clone("FDAcceptanceCorrection_True_HEURSCALE"));
  TH1D *LinCombOneShotETrue_POTScaled_HEURSCALE = static_cast<TH1D *>(
      LinCombOneShotETrue->Clone("LinCombOneShotETrue_POTScaled_HEURSCALE"));

  FDToNDScale = LinCombOneShotETrue_POTScaled_HEURSCALE->GetMaximum() /
                       FDPrediction_True_POTScaled_HEURSCALE->GetMaximum();

  FDPrediction_True_POTScaled_HEURSCALE->Scale(FDToNDScale);
  FDFluxCorrection_True_HEURSCALE->Scale(FDToNDScale);
  FDAcceptanceCorrection_True_HEURSCALE->Scale(FDToNDScale);

  TH1D *LinCombOneShotETrue_POTScaled_WithCorr_HEURSCALE =
      static_cast<TH1D *>(LinCombOneShotETrue_POTScaled_HEURSCALE->Clone(
          "LinCombOneShotETrue_POTScaled_WithCorr_HEURSCALE"));

  LinCombOneShotETrue_POTScaled_WithCorr_HEURSCALE->Add(
      FDFluxCorrection_True_HEURSCALE);
  LinCombOneShotETrue_POTScaled_WithCorr_HEURSCALE->Add(
      FDAcceptanceCorrection_True_HEURSCALE);

  FDPrediction_True_POTScaled->SetDirectory(resD);
  FDFluxCorrection_True_clone->SetDirectory(resD);
  FDAcceptanceCorrection_True_clone->SetDirectory(resD);
  LinCombOneShotETrue_POTScaled->SetDirectory(resD);

  outfile->Write();
}
