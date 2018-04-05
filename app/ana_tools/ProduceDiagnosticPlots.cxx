#include "DepositsSummaryTreeReader.hxx"
#include "StopConfigTreeReader.hxx"
#include "SliceConfigTreeReader.hxx"
#include "SelectionSummaryTreeReader.hxx"
#include "SimConfigTreeReader.hxx"

#include "EffCorrector.hxx"

#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

#include "BoundingBox.hxx"

#include "TTree.h"

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>

std::string InputDepostisSummaryFile;
std::string OutputFile;
std::string EffCorrectorFile;
std::string FluxFitFile;

std::string DefaultEnergyBinning = "0_10:0.1";

double HadrVeto = 0;  // GeV
bool SelMuExit = false;
double SelMuExitKE = 0;  // GeV
double HadrVisAcceptance = 0; //GeV
std::vector<double> ExtraVertexSelectionPadding;

enum SelLevel { kUnSelected = 0, kSelectedMu,
  kSelectedHadr, kSelected, kCorrected, kNSelLevels };

std::string to_str(SelLevel sl){
  switch(sl){
    case kUnSelected:{
      return "UnSelected";
    }
    case kSelectedMu:{
      return "SelectedMu";
    }
    case kSelectedHadr:{
      return "SelectedHadr";
    }
    case kSelected:{
      return "Selected";
    }
    case kCorrected:{
      return "Corrected";
    }
    default: { throw; }
  }
}

SelLevel GetSelectedLevel(DepositsSummary &edr, std::vector<BoundingBox> &FVs){

  bool NumuInStop = ((edr.stop >= 0) && (edr.PrimaryLepPDG == 13));
  bool InFV = FVs[edr.stop].Contains({edr.vtx[0], edr.vtx[1], edr.vtx[2]});
  if(!(NumuInStop&&InFV)){
    return kUnSelected;
  }

  bool SelMu = ((SelMuExitKE == 0) || (edr.LepExitKE > SelMuExitKE));
  bool SelHadr = (edr.TotalNonlep_Dep_veto < HadrVeto);
  bool HadrAccepted = ((HadrVisAcceptance == 0) ||
    (edr.GetProjection(DepositsSummary::kEHadr_vis) <= HadrVisAcceptance));

  if(SelMu && SelHadr && HadrAccepted){
    return kSelected;
  }

  if(SelHadr && HadrAccepted){
    return kSelectedHadr;
  } else {
    return kSelectedMu;
  }
}

size_t gi(SelLevel sl){
  return static_cast<size_t>(sl);
}
SelLevel gsl(size_t si){
  return static_cast<SelLevel>(si);
}

EffCorrector::ModeEnum EffCorrMode = EffCorrector::kEHadrVisDetPos;

std::vector<double> GPVB(DepositsSummary::ProjectionVar pv){
  static std::map<DepositsSummary::ProjectionVar, std::vector<double> > cache;

  if(cache.count(pv)){
    return cache[pv];
  }

  switch(pv){
    case DepositsSummary::kETrue:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kEAvail_True:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kEHadr_True:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kEFSLep_True:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kEHadr_vis:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kEHadrLate_vis:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kEHadrAll_vis:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kEFSLep_vis:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kEFSLepLate_vis:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kEFSLepAll_vis:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kEFSLepAndDescendents_vis:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kEFSLepAndDescendentsLate_vis:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kEFSLepAndDescendentsAll_vis:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kERec:{
      cache[pv] = BuildBinEdges(DefaultEnergyBinning);
      break;
    }
    case DepositsSummary::kERecResidual:{
      cache[pv] = BuildBinEdges("-2_2:0.025");
      break;
    }
    case DepositsSummary::kERecBias:{
      cache[pv] = BuildBinEdges("-1_1:0.02");
      break;
    }
    case DepositsSummary::kNProjectionVars:
    default: {
      cache[pv] = std::vector<double>{};
    }
  }
  return GPVB(pv);
}

void SayUsage(char const *argv[]) {
  std::cout
      << "[USAGE]: " << argv[0]
      << "\n"
         "\t-i <stopprocessor.root>     : TChain descriptor for"
         " input tree. \n"
         "\t-o <outputfile.root>        : Output file to write "
         "selected tree to.\n"
         "\t-v <hadr veto threshold>    : Hadronic shower veto threshold in "
         "MeV.\n"
         "\t-FV <fvx,y,z>               : Vertex selection fiducial volume "
         "padding inside of the \n"
         "\t                              non-veto active region "
         "{Default: 0,0,0}.\n"
         "\t-m <muon exit KE>           : Muon exit threshold KE in MeV.\n"
         "\t-A <EHadrVisMax_GeV>        : Acceptance cut for hadronic shower"
         " energy. Showers with\n"
         "\t                              more than this energy are cut and "
         "filled in by far detector MC.\n"
         "\t-E <effcorrector.root>      : File containing efficiency "
                  "histograms.\n"
         "\t-F <fluxfitresults.root>        : Result of "
          "dp_FitFluxes to use to build observation.\n"
          "\t-M <eff correction mode>    : Kinematics to use to determine the "
          "hadronic shower selection\n"
          "\t                              efficiency correction. {Default = 3}"
          "\n"
          "\t                              1: True hadronic energy, absolute "
          "off axis position.\n"
          "\t                              2: True hadronic energy, position "
          "within the detector.\n"
          "\t                              3: Visible hadronic energy, position "
          "within the detector.\n"
          "\t-b <binning descriptor>     : Energy binning descriptor that can "
          "be used to produce a perfect efficiency correction in the relevant "
          "off-axis bins.\n"
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
      InputDepostisSummaryFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-o") {
      OutputFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-v") {
      HadrVeto = str2T<double>(argv[++opt]) * 1E-3;
    } else if (std::string(argv[opt]) == "-m") {
      SelMuExitKE = str2T<double>(argv[++opt]) * 1E-3;
      SelMuExit = true;
    } else if (std::string(argv[opt]) == "-A") {
      HadrVisAcceptance = str2T<double>(argv[++opt]);
    }
    else if (std::string(argv[opt]) == "-FV") {
      ExtraVertexSelectionPadding = ParseToVect<double>(argv[++opt], ",");
      if(ExtraVertexSelectionPadding.size() != 3){
        std::cout << "[ERROR]: -FV option contained " <<
          ExtraVertexSelectionPadding.size() << " entries, expected 3."
          << std::endl;
        SayUsage(argv);
        exit(1);
      }
    } else if (std::string(argv[opt]) == "-E") {
      EffCorrectorFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-F") {
      FluxFitFile = argv[++opt];
    } else if (std::string(argv[opt]) == "-M") {
      EffCorrMode = static_cast<EffCorrector::ModeEnum>(str2T<int>(argv[++opt]));
    } else if (std::string(argv[opt]) == "-b") {
      DefaultEnergyBinning = argv[++opt];
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

  if(!ExtraVertexSelectionPadding.size()){
      int argc_dum = 3;
      char const * argv_dum[] = { argv[0], "-FV", "0,0,0"};
      handleOpts(argc_dum, argv_dum);
  }

  SimConfig simCRdr(InputDepostisSummaryFile);
  StopConfig csRdr(InputDepostisSummaryFile);

  std::vector<BoundingBox> FVs = csRdr.GetStopBoundingBoxes(true,
    {ExtraVertexSelectionPadding[0], ExtraVertexSelectionPadding[1],
      ExtraVertexSelectionPadding[2]});

  SliceConfig slCfg(FluxFitFile);
  std::vector< std::pair<double, double> > XRange_comp = slCfg.GetXRanges();
  std::vector<double> Coeffs_comp = slCfg.GetCoeffs();

  std::pair< std::vector<double>, std::vector<double> >
    xrb = SliceConfig::BuildXRangeBinsCoeffs(XRange_comp, Coeffs_comp.data(),
      true);

  std::vector<double> &XRangeBins = xrb.first;

  TH1D *SliceHelper = new TH1D(
      "SliceHelper", "", (XRangeBins.size() - 1), XRangeBins.data());

  TFile *of = CheckOpenFile(OutputFile, "RECREATE");

  //Stop Plots
  std::array<
    std::map<
      DepositsSummary::ProjectionVar,
      std::vector<TH1D *> >, kNSelLevels> Stop_plots;
  std::array<std::vector<TH2D *>, kNSelLevels> Stop_ERecETrue;
  //Slice Plots
  std::array<
    std::map<
      DepositsSummary::ProjectionVar,
      std::vector<TH1D *> >, kNSelLevels> Slice_plots;
  std::array<std::vector<TH2D *>, kNSelLevels> Slice_ERecETrue;
  //integrated
  std::array<
    std::map<
      DepositsSummary::ProjectionVar, TH1D * >, kNSelLevels> Integrated_plots;
  std::array<TH2D *, kNSelLevels> Integrated_ERecETrue;
  //AbsPos
  std::array<
    std::map<
      DepositsSummary::ProjectionVar, TH2D * >, kNSelLevels> AbsPos_plots;
  std::array<TH1D *, kNSelLevels> AbsPos_plots1D;

  TDirectory *StopDir = of->mkdir("StopPlots");
  TDirectory *SliceDir = of->mkdir("SlicePlots");
  TDirectory *IntegratedDir = of->mkdir("IntegratedPlots");

  for(size_t sel_it = gi(kUnSelected); sel_it < gi(kNSelLevels); ++sel_it) {
    //mv to dir
    StopDir->cd();
    for(size_t stop_it = 0; FVs.size() && (stop_it < FVs.size()); ++stop_it){
      DepositsSummary::ProjectionVar pv = DepositsSummary::kETrue;
      while(pv != DepositsSummary::kNProjectionVars){
        std::vector<double> const & binning = GPVB(pv);
        if(!binning.size()){
          std::cout << "Failed to get binning for " << to_str(pv) << std::endl;
          throw;
        }

        Stop_plots[sel_it][pv].push_back(new TH1D(
          (std::string("Stop_") + to_str(stop_it) + "_" + to_str(pv) + "_"
            + to_str(gsl(sel_it))).c_str(),
          (std::string(";") + to_title(pv) +
            ";Count").c_str(),
          (binning.size() - 1), binning.data() ));
        pv = next(pv);
      }
      Stop_ERecETrue[sel_it].push_back(
        new TH2D((std::string("Stop_") + to_str(stop_it) + "_"
          + to_str(gsl(sel_it)) + "_ERecETrue").c_str(),
          (std::string(";") + to_title(DepositsSummary::kETrue)
          + to_title(DepositsSummary::kERec) + ";Count").c_str(),
        (GPVB(DepositsSummary::kETrue).size() - 1),
          GPVB(DepositsSummary::kETrue).data(),
        (GPVB(DepositsSummary::kERec).size() - 1),
          GPVB(DepositsSummary::kERec).data() ));
    }
    SliceDir->cd();
    for(size_t slice_it = 0; XRangeBins.size() && (slice_it < XRangeBins.size());
      ++slice_it){

      DepositsSummary::ProjectionVar pv = DepositsSummary::kETrue;
      while(pv != DepositsSummary::kNProjectionVars){
        Slice_plots[sel_it][pv].push_back(new TH1D(
          (std::string("Slice_") + to_str(slice_it) + "_" + to_str(pv) + "_"
            + to_str(gsl(sel_it))).c_str(),
          (std::string(";") + to_title(pv) +
            ";Count").c_str(),
          (GPVB(pv).size() - 1), GPVB(pv).data() ));
        pv = next(pv);
      }

      Slice_ERecETrue[sel_it].push_back(
        new TH2D((std::string("Slice_") + to_str(slice_it) + "_"
          + to_str(gsl(sel_it)) + "_ERecETrue").c_str(),
        (std::string(";") + to_title(DepositsSummary::kETrue)
          + to_title(DepositsSummary::kERec) + ";Count").c_str(),
          (GPVB(DepositsSummary::kETrue).size() - 1),
            GPVB(DepositsSummary::kETrue).data(),
          (GPVB(DepositsSummary::kERec).size() - 1),
            GPVB(DepositsSummary::kERec).data() ));

    }
    IntegratedDir->cd();
    DepositsSummary::ProjectionVar pv = DepositsSummary::kETrue;
    while(pv != DepositsSummary::kNProjectionVars){
      Integrated_plots[sel_it][pv] = new TH1D(
        ("Integrated_" + to_str(pv) + "_"
          + to_str(gsl(sel_it))).c_str(),
        (std::string(";") + to_title(pv) +
          ";Count").c_str(),
        (GPVB(pv).size() - 1), GPVB(pv).data() );

      AbsPos_plots[sel_it][pv] = new TH2D((
        std::string("AbsPos_") + to_str(pv) + "_"
          + to_str(gsl(sel_it))).c_str(),
          (std::string(";") + to_title(pv)
          + "Off-axis position (cm);Count").c_str(),
          (GPVB(pv).size() - 1), GPVB(pv).data(),
          400, -3800, 200 );
      pv = next(pv);
    }
    Integrated_ERecETrue[sel_it] = new TH2D(
      (std::string("Integrated_ERecETrue_") + to_str(gsl(sel_it))).c_str(),
    (std::string(";") + to_title(DepositsSummary::kETrue)
      + to_title(DepositsSummary::kERec) + ";Count").c_str(),
      (GPVB(DepositsSummary::kETrue).size() - 1),
        GPVB(DepositsSummary::kETrue).data(),
      (GPVB(DepositsSummary::kERec).size() - 1),
        GPVB(DepositsSummary::kERec).data() );
    AbsPos_plots1D[sel_it]  = new TH1D((
        std::string("AbsPos_") + "_" + to_str(gsl(sel_it))).c_str(),
        (std::string(";") + "Off-axis position (cm);Count").c_str(),
        400, -3800, 200 );

  }

  EffCorrector EffCorr(EffCorrMode, EffCorrectorFile, csRdr);

  DepositsSummary DepSumRdr(InputDepostisSummaryFile);
  size_t loud_every = DepSumRdr.GetEntries() / 10;
  Long64_t NEntries = DepSumRdr.GetEntries();
  for (Long64_t e_it = 0; e_it < NEntries; ++e_it) {
    DepSumRdr.GetEntry(e_it);
    if (loud_every && !(e_it % loud_every)) {
      std::cout << "\r[INFO]: Read " << e_it << " entries... ( vtx: {"
                << DepSumRdr.vtx[0] << ", " << DepSumRdr.vtx[1] << ", "
                << DepSumRdr.vtx[2] << "}, Enu: " << DepSumRdr.nu_4mom[3]
                << " )" << std::endl;
    }

    SelLevel sl = GetSelectedLevel(DepSumRdr, FVs);

    Int_t slice =
        SliceHelper->GetXaxis()->FindFixBin(DepSumRdr.vtx[0]) - 1;
    if(slice >= Int_t(XRangeBins.size())){
      slice = -1;
    }

    double mueffweight = 1;
    if (SelMuExit) {
      mueffweight = EffCorr.GetMuonKinematicsEffWeight(DepSumRdr);
    }
    double hadreffweight = EffCorr.GetHadronKinematicsEffWeight(DepSumRdr);
    double effweight = mueffweight * hadreffweight;

    DepositsSummary::ProjectionVar pv = DepositsSummary::kETrue;
    while(pv != DepositsSummary::kNProjectionVars){
      Stop_plots[gi(kUnSelected)][pv][DepSumRdr.stop]->Fill(
        DepSumRdr.GetProjection(pv), DepSumRdr.stop_weight);
      if(slice != -1){
        Slice_plots[gi(kUnSelected)][pv][slice]->Fill(DepSumRdr.GetProjection(pv),
          DepSumRdr.stop_weight);
      }
      Integrated_plots[gi(kUnSelected)][pv]->Fill(DepSumRdr.GetProjection(pv),
        DepSumRdr.stop_weight);
      AbsPos_plots[gi(kUnSelected)][pv]->Fill(DepSumRdr.GetProjection(pv),
        DepSumRdr.vtx[0], DepSumRdr.stop_weight);

      if(sl == kUnSelected){
        pv = next(pv);
        continue;
      }

      Stop_plots[gi(sl)][pv][DepSumRdr.stop]->Fill(
        DepSumRdr.GetProjection(pv), DepSumRdr.stop_weight);
      if(slice != -1){
        Slice_plots[gi(sl)][pv][slice]->Fill(DepSumRdr.GetProjection(pv),
          DepSumRdr.stop_weight);
      }
      Integrated_plots[gi(sl)][pv]->Fill(DepSumRdr.GetProjection(pv),
        DepSumRdr.stop_weight);
      AbsPos_plots[gi(sl)][pv]->Fill(DepSumRdr.GetProjection(pv),
        DepSumRdr.vtx[0], DepSumRdr.stop_weight);

      if(sl != kSelected){
        pv = next(pv);
        continue;
      }

      Stop_plots[gi(kSelectedMu)][pv][DepSumRdr.stop]->Fill(
        DepSumRdr.GetProjection(pv), DepSumRdr.stop_weight);
      if(slice != -1){
        Slice_plots[gi(kSelectedMu)][pv][slice]->Fill(DepSumRdr.GetProjection(pv),
          DepSumRdr.stop_weight);
      }
      Integrated_plots[gi(kSelectedMu)][pv]->Fill(DepSumRdr.GetProjection(pv),
        DepSumRdr.stop_weight);
      AbsPos_plots[gi(kSelectedMu)][pv]->Fill(DepSumRdr.GetProjection(pv),
        DepSumRdr.vtx[0], DepSumRdr.stop_weight);

      Stop_plots[gi(kSelectedHadr)][pv][DepSumRdr.stop]->Fill(
        DepSumRdr.GetProjection(pv), DepSumRdr.stop_weight);
      if(slice != -1){
        Slice_plots[gi(kSelectedHadr)][pv][slice]->Fill(DepSumRdr.GetProjection(pv),
          DepSumRdr.stop_weight);
      }
      Integrated_plots[gi(kSelectedHadr)][pv]->Fill(DepSumRdr.GetProjection(pv),
        DepSumRdr.stop_weight);
      AbsPos_plots[gi(kSelectedHadr)][pv]->Fill(DepSumRdr.GetProjection(pv),
        DepSumRdr.vtx[0], DepSumRdr.stop_weight);


      Stop_plots[gi(kCorrected)][pv][DepSumRdr.stop]->Fill(
        DepSumRdr.GetProjection(pv), DepSumRdr.stop_weight * effweight);
      if(slice != -1){
        Slice_plots[gi(kCorrected)][pv][slice]->Fill(DepSumRdr.GetProjection(pv),
          DepSumRdr.stop_weight * effweight);
      }
      Integrated_plots[gi(kCorrected)][pv]->Fill(DepSumRdr.GetProjection(pv),
        DepSumRdr.stop_weight * effweight);
      AbsPos_plots[gi(kCorrected)][pv]->Fill(DepSumRdr.GetProjection(pv),
        DepSumRdr.vtx[0], DepSumRdr.stop_weight * effweight);

      pv = next(pv);
    }

    Stop_ERecETrue[gi(kUnSelected)][DepSumRdr.stop]->Fill(
      DepSumRdr.GetProjection(DepositsSummary::kETrue),
      DepSumRdr.GetProjection(DepositsSummary::kERec),
      DepSumRdr.stop_weight);
    if(slice != -1){
      Slice_ERecETrue[gi(kUnSelected)][slice]->Fill(
        DepSumRdr.GetProjection(DepositsSummary::kETrue),
        DepSumRdr.GetProjection(DepositsSummary::kERec),
        DepSumRdr.stop_weight);
    }
    Integrated_ERecETrue[gi(kUnSelected)]->Fill(
      DepSumRdr.GetProjection(DepositsSummary::kETrue),
      DepSumRdr.GetProjection(DepositsSummary::kERec),
      DepSumRdr.stop_weight);
    AbsPos_plots1D[gi(kUnSelected)]->Fill(DepSumRdr.vtx[0],
      DepSumRdr.stop_weight);

    if(sl == kUnSelected){
      continue;
    }

    Stop_ERecETrue[gi(sl)][DepSumRdr.stop]->Fill(
      DepSumRdr.GetProjection(DepositsSummary::kETrue),
      DepSumRdr.GetProjection(DepositsSummary::kERec),
      DepSumRdr.stop_weight);
    if(slice != -1){
      Slice_ERecETrue[gi(sl)][slice]->Fill(
        DepSumRdr.GetProjection(DepositsSummary::kETrue),
        DepSumRdr.GetProjection(DepositsSummary::kERec),
        DepSumRdr.stop_weight);
    }
    Integrated_ERecETrue[gi(sl)]->Fill(
      DepSumRdr.GetProjection(DepositsSummary::kETrue),
      DepSumRdr.GetProjection(DepositsSummary::kERec),
      DepSumRdr.stop_weight);
    AbsPos_plots1D[gi(sl)]->Fill(DepSumRdr.vtx[0], 
      DepSumRdr.stop_weight);

    if(sl != kSelected){
      continue;
    }

    Stop_ERecETrue[gi(kSelectedMu)][DepSumRdr.stop]->Fill(
      DepSumRdr.GetProjection(DepositsSummary::kETrue),
      DepSumRdr.GetProjection(DepositsSummary::kERec),
      DepSumRdr.stop_weight);
    if(slice != -1){
      Slice_ERecETrue[gi(kSelectedMu)][slice]->Fill(
        DepSumRdr.GetProjection(DepositsSummary::kETrue),
        DepSumRdr.GetProjection(DepositsSummary::kERec),
        DepSumRdr.stop_weight);
    }
    Integrated_ERecETrue[gi(kSelectedMu)]->Fill(
      DepSumRdr.GetProjection(DepositsSummary::kETrue),
      DepSumRdr.GetProjection(DepositsSummary::kERec),
      DepSumRdr.stop_weight);
    AbsPos_plots1D[gi(kSelectedMu)]->Fill(DepSumRdr.vtx[0],
      DepSumRdr.stop_weight);

    Stop_ERecETrue[gi(kSelectedHadr)][DepSumRdr.stop]->Fill(
      DepSumRdr.GetProjection(DepositsSummary::kETrue),
      DepSumRdr.GetProjection(DepositsSummary::kERec),
      DepSumRdr.stop_weight);
    if(slice != -1){
      Slice_ERecETrue[gi(kSelectedHadr)][slice]->Fill(
        DepSumRdr.GetProjection(DepositsSummary::kETrue),
        DepSumRdr.GetProjection(DepositsSummary::kERec),
        DepSumRdr.stop_weight);
    }
    Integrated_ERecETrue[gi(kSelectedHadr)]->Fill(
      DepSumRdr.GetProjection(DepositsSummary::kETrue),
      DepSumRdr.GetProjection(DepositsSummary::kERec),
      DepSumRdr.stop_weight);
    AbsPos_plots1D[gi(kSelectedHadr)]->Fill(DepSumRdr.vtx[0],
      DepSumRdr.stop_weight);

    Stop_ERecETrue[gi(kCorrected)][DepSumRdr.stop]->Fill(
      DepSumRdr.GetProjection(DepositsSummary::kETrue),
      DepSumRdr.GetProjection(DepositsSummary::kERec),
      DepSumRdr.stop_weight * effweight);
    if(slice != -1){
      Slice_ERecETrue[gi(kCorrected)][slice]->Fill(
        DepSumRdr.GetProjection(DepositsSummary::kETrue),
        DepSumRdr.GetProjection(DepositsSummary::kERec),
        DepSumRdr.stop_weight * effweight);
    }
    Integrated_ERecETrue[gi(kCorrected)]->Fill(
      DepSumRdr.GetProjection(DepositsSummary::kETrue),
      DepSumRdr.GetProjection(DepositsSummary::kERec),
      DepSumRdr.stop_weight * effweight);
    AbsPos_plots1D[gi(kCorrected)]->Fill(DepSumRdr.vtx[0],
      DepSumRdr.stop_weight * effweight);
  }

  of->Write();
}
