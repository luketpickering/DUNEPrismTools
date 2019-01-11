#define FLS_WRAP_IN_NAMESPACE
#include "FluxLinearSolver_Standalone.hxx"

#include "CAFReader.hxx"
#include "PlotProcessor.hxx"
#include "ROOTUtility.hxx"

#include "DUNETDRDetHelper.hxx"
#include "OscillationHelper.hxx"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

struct OAWeigher {
  Eigen::VectorXd OACoeffs;
  double FVRatio;

  double GetWeight(double x) {
    Int_t bin = std::floor(x * 2.0);
    if (bin < 0) {
      return 0;
    }
    if (bin >= OACoeffs.size()) {
      return 0;
    }
    // std::cout << "[INFO]: x = " << x << " m, gives bin " << bin
    //           << ", FVRatio = " << FVRatio << std::endl;
    return OACoeffs[bin] * FVRatio;
  }
};

// Expect a fhicl like
// Input_ND: ND_CAF.root
// Input_FD: FD_CAF.root
// Oscillation:
// Input_ND_Flux: { File: name.root Hist: name }
// Input_FD_Flux: { File: name.root Hist: name }
// ETrueBinning: [<nbins>,<histmin>,<histmax>]
int main(int argc, char const *argv[]) {

  TH1::SetDefaultSumw2();

  if (argc != 2) {
    std::cout << "[ERROR]: Expected to be passed a single FHiCL file "
                 "describing the inputs."
              << std::endl;
    return 1;
  }

  fhicl::ParameterSet ps = fhicl::make_ParameterSet(argv[1]);

  std::array<double, 3> ETrueBinning =
      ps.get<std::array<double, 3>>("ETrueBinning");

  fls::FluxLinearSolver linsolv;
  fls::FluxLinearSolver::Params solver_params =
      fls::FluxLinearSolver::GetDefaultParams();
  solver_params.algo_id = fls::FluxLinearSolver::Params::kInverse;
  solver_params.OORMode = fls::FluxLinearSolver::Params::kGaussianDecay;
  solver_params.OORSide = fls::FluxLinearSolver::Params::kLeft;
  solver_params.ExpDecayRate = 3;
  solver_params.OORFactor = 0.1;
  solver_params.FitBetweenFoundPeaks = true;
  solver_params.MergeENuBins = 2;
  solver_params.OffAxisRangesDescriptor = "0_32:0.5";

  std::string NDFluxFile =
      ps.get<fhicl::ParameterSet>("Input_ND_Flux").get<std::string>("File");
  std::string NDFluxHist =
      ps.get<fhicl::ParameterSet>("Input_ND_Flux").get<std::string>("Hist");
  std::string FDFluxFile =
      ps.get<fhicl::ParameterSet>("Input_FD_Flux").get<std::string>("File");
  std::string FDFluxHist =
      ps.get<fhicl::ParameterSet>("Input_FD_Flux").get<std::string>("Hist");

  linsolv.Initialize(solver_params, {NDFluxFile, NDFluxHist},
                     {FDFluxFile, FDFluxHist});

  OscillationHelper oh;
  oh.Setup(ps.get<fhicl::ParameterSet>("Oscillation"));

  std::array<double, 6> oscparams;
  std::copy_n(oh.OscParams, 6, oscparams.begin());
  linsolv.OscillateFDFlux(oscparams, {14, 14}, oh.DipAngle_degrees);
  Eigen::VectorXd cv = linsolv.Solve(1E-9);
  OAWeigher weighter;

  weighter.OACoeffs = cv;
  weighter.FVRatio = FD_ND_FVRatio(50);

  CAFReader ND_rdr(ps.get<std::string>("Input_ND"));
  CAFReader FD_rdr(ps.get<std::string>("Input_FD"));

  Plot1DCAF FD_ETrue(
      [](CAFReader const &ev) -> std::array<double, 2> {
        return {ev.Ev, ev.POTWeight};
      },
      "Enu_true", ";E_{#nu} (GeV); Count/POT", ETrueBinning[0], ETrueBinning[1],
      ETrueBinning[2]);

  size_t FDEvents = FD_rdr.GetEntries();
  for (size_t fd_it = 0; fd_it < FDEvents; ++fd_it) {
    FD_rdr.GetEntry(fd_it);
    if (!FV_Select(FD_rdr)) {
      continue;
    }

    oh.SetOscillationChannel(FD_rdr.nuPDGunosc, FD_rdr.nuPDG);
    double ow = oh.GetWeight(FD_rdr.Ev);

    bool true_sel = FHC_Numu_True_Select(FD_rdr);
    if (true_sel) {
      FD_ETrue.Process(FD_rdr, ow);
    }

    // if (FD_FHC_Numu_Select(FD_rdr)) {
    // }
  }

  Plot2DCAF NDEventRate(
      [](CAFReader const &ev) -> std::array<double, 3> {
        return {ev.Ev, ev.det_x + ev.vtx_x * 1E-2, ev.POTWeight};
      },
      "OffAxisEvRateEnu", "", ETrueBinning[0], ETrueBinning[1], ETrueBinning[2],
      90, -5.25, 39.75);

  Plot2DCAF NDEventRate_fluxbin(
      [](CAFReader const &ev) -> std::array<double, 3> {
        return {ev.Ev, ev.det_x + ev.vtx_x * 1E-2, ev.POTWeight};
      },
      "OffAxisEvRateEnu_fluxbin", "", 400, 0, 10, 91, -0.25, 45.25);

  Plot2DCAF NDEventRate_lcw(NDEventRate.Clone("OffAxisEvRateEnu_weighted"));

  size_t NDEvents = ND_rdr.GetEntries();
  for (size_t nd_it = 0; nd_it < NDEvents; ++nd_it) {
    ND_rdr.GetEntry(nd_it);
    if (!FV_Select(ND_rdr)) {
      continue;
    }

    if (FHC_Numu_True_Select(ND_rdr)) {
      NDEventRate.Process(ND_rdr);
      NDEventRate_fluxbin.Process(ND_rdr);

      NDEventRate_lcw.Process(
          ND_rdr, weighter.GetWeight(ND_rdr.det_x + ND_rdr.vtx_x * 1.E-2));
    }
    // if (ND_FHC_Numu_Select(ND_rdr)) {
    // }
  }

  TH2D *NDFlux = GetHistogram<TH2D>(NDFluxFile, NDFluxHist);

  for (Int_t ybin = 0; ybin < NDFlux->GetYaxis()->GetNbins(); ++ybin) {
    double w = weighter.GetWeight(NDFlux->GetYaxis()->GetBinCenter(ybin + 1));
    for (Int_t xbin = 0; xbin < NDFlux->GetXaxis()->GetNbins(); ++xbin) {
      NDFlux->SetBinContent(xbin + 1, ybin + 1,
                            NDFlux->GetBinContent(xbin + 1, ybin + 1) * w);
    }
  }

  std::string oupf = ps.get<std::string>("OutputFile");
  TFile *oupFile = CheckOpenFile(oupf, "RECREATE");

  TDirectory *fit_dir = oupFile->mkdir("FluxFit");

  linsolv.Write(fit_dir);

  NDFlux->SetDirectory(fit_dir);
  TH1D *NDF_prjx = NDFlux->ProjectionX("FluxLinComb");
  NDF_prjx->SetDirectory(fit_dir);

  TDirectory *evr_dir = oupFile->mkdir("EventRates");

  FD_ETrue.Write(evr_dir);
  NDEventRate.Write(evr_dir);
  NDEventRate_fluxbin.Write(evr_dir);

  TH1D *prjx = NDEventRate_lcw.hist->ProjectionX("LinComb");
  prjx->SetDirectory(evr_dir);
  NDEventRate_lcw.Write(evr_dir);

  oupFile->Write();
  oupFile->Close();
}
