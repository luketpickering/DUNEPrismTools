#define FLS_WRAP_IN_NAMESPACE
#include "FluxLinearSolver_Standalone.hxx"

#include "CovarianceHelper.hxx"
#include "EffectiveFluxUncertaintyHelper.hxx"
#include "OscillationHelper.hxx"
#include "ROOTUtility.hxx"

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"

#include <memory>
#include <random>

// Expect a fhicl like
// Oscillation: {}
// Input_ND_Flux: { File: name.root Hist: name }
// Input_FD_Flux: { File: name.root Hist: name }
// Uncertainties: {}
int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();

  if (argc != 2) {
    std::cout << "[ERROR]: Expected to be passed a single FHiCL file "
                 "describing the inputs."
              << std::endl;
    return 1;
  }
  fhicl::ParameterSet ps = fhicl::make_ParameterSet(argv[1]);

  fls::FluxLinearSolver linsolv;
  fls::FluxLinearSolver::Params solver_params =
      fls::FluxLinearSolver::GetDefaultParams();
  solver_params.algo_id = fls::FluxLinearSolver::Params::kInverse;
  solver_params.OORMode = fls::FluxLinearSolver::Params::kGaussianDecay;
  solver_params.OORSide = fls::FluxLinearSolver::Params::kLeft;
  solver_params.ExpDecayRate = 3;
  solver_params.OORFactor = 0.1;
  solver_params.FitBetweenFoundPeaks = false;
  solver_params.FitBetween =
      ps.get<std::pair<double, double>>("FitBetween", {0.5, 4});
  solver_params.MergeENuBins = 0;
  solver_params.OffAxisRangesDescriptor = "0_57:1";

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

  size_t NThrows = ps.get<size_t>("NThrows");

  EffectiveFluxUncertaintyHelper flux_uncert;
  flux_uncert.Initialize(ps.get<fhicl::ParameterSet>("Uncertainties"));

  size_t NFluxParameters = flux_uncert.GetNParameters();
  std::vector<double> flux_param_values(NFluxParameters);
  std::random_device r;

  std::mt19937 RNEngine(r());
  std::normal_distribution<double> RNJesus(0, 1);
  std::unique_ptr<TH2> NDFluxes_nom(
      dynamic_cast<TH2 *>(linsolv.NDFluxes->Clone()));
  NDFluxes_nom->SetDirectory(nullptr);
  std::unique_ptr<TH2> NDFluxes_var(
      dynamic_cast<TH2 *>(linsolv.NDFluxes->Clone()));
  NDFluxes_var->SetDirectory(nullptr);

  std::unique_ptr<TH1> FDFlux_nom(
      dynamic_cast<TH1 *>(linsolv.FDFlux_unosc->Clone()));
  FDFlux_nom->SetDirectory(nullptr);
  std::unique_ptr<TH1> FDFlux_var(
      dynamic_cast<TH1 *>(linsolv.FDFlux_unosc->Clone()));
  FDFlux_var->SetDirectory(nullptr);

  std::string oupf = ps.get<std::string>("OutputFile");
  TFile *oupFile = CheckOpenFile(oupf, "RECREATE");

  bool SaveThrows = ps.get<bool>("SaveThrows", false);

  int nu_config_nd = flux_uncert.GetNuConfig(14, true, true);
  int nu_config_fd = flux_uncert.GetNuConfig(14, false, true);

  std::unique_ptr<TH2> ThrowCoefficients(
      new TH2D("ThrowCoefficients", ";Off axis angle (mrads);Coefficient", 60,
               0, 60, 1000, -1E-7, 1E-7));
  std::vector<std::vector<double>> RelativeFluxMatches;
  std::unique_ptr<TH2> RelativeFluxDifferences(new TH2D(
      "RelativeFluxDifferences", ";E_{#nu} (GeV);(LinComb - FDOsc)/FDUnosc",
      FDFlux_nom->GetXaxis()->GetNbins(),
      FDFlux_nom->GetXaxis()->GetBinLowEdge(1),
      FDFlux_nom->GetXaxis()->GetBinUpEdge(FDFlux_nom->GetXaxis()->GetNbins()),
      100, -2, 2));
  RelativeFluxDifferences->SetDirectory(nullptr);

  for (size_t t_it = 0; t_it < NThrows; ++t_it) {

    if (t_it && !(t_it % 100)) {
      std::cout << "[INFO]: Done " << t_it << " throws." << std::endl;
    }

    // Throw flux parameters
    for (double &pval : flux_param_values) {
      pval = RNJesus(RNEngine);
    }

    for (Int_t enu_bin_it = 0;
         enu_bin_it < NDFluxes_nom->GetXaxis()->GetNbins(); ++enu_bin_it) {
      double enu_GeV = NDFluxes_nom->GetXaxis()->GetBinCenter(enu_bin_it + 1);
      for (Int_t OA_bin_it = 0;
           OA_bin_it < NDFluxes_nom->GetYaxis()->GetNbins(); ++OA_bin_it) {
        double oa_mrad = NDFluxes_nom->GetYaxis()->GetBinCenter(OA_bin_it + 1);

        // Reweight NDfluxes
        {
          int bin_nd = flux_uncert.GetBin(14, enu_GeV, oa_mrad, 0, true);
          double weight = 1;
          size_t p_it = 0;
          for (p_it = 0; p_it < NFluxParameters; ++p_it) {
            weight *= flux_uncert.GetFluxWeight(p_it, flux_param_values[p_it],
                                                bin_nd, nu_config_nd);
          }

          NDFluxes_var->SetBinContent(
              enu_bin_it + 1, OA_bin_it + 1,
              NDFluxes_nom->GetBinContent(enu_bin_it + 1, OA_bin_it + 1) *
                  weight);
        }
      }

      // ReWeight FDFlux
      {

        int bin_fd = flux_uncert.GetBin(14, enu_GeV, 0, 0, false);

        double weight = 1;
        size_t p_it = 0;
        for (p_it = 0; p_it < NFluxParameters; ++p_it) {
          weight *= flux_uncert.GetFluxWeight(p_it, flux_param_values[p_it],
                                              bin_fd, nu_config_fd);
        }

        FDFlux_var->SetBinContent(
            enu_bin_it + 1,
            FDFlux_nom->GetBinContent(enu_bin_it + 1) * weight);
      }
    }

    linsolv.SetNDFluxes(NDFluxes_var.get());
    linsolv.SetFDFluxUnOsc(FDFlux_var.get());
    linsolv.OscillateFDFlux(oscparams, {14, 14}, oh.DipAngle_degrees);
    Eigen::VectorXd cv = linsolv.Solve(1E-9);
    // Add to TH2 of coeffs and best fit
    for (int it = 0; it < cv.size(); ++it) {
      ThrowCoefficients->Fill(double(it) + 0.5, cv[it]);
    }

    // Add to best fit covmat.
    std::unique_ptr<TH1> Match(linsolv.GetLastMatch());
    std::unique_ptr<TH1> FDFlux_unosc(
        dynamic_cast<TH1 *>(linsolv.FDFlux_unosc->Clone()));
    FDFlux_unosc->SetDirectory(nullptr);
    std::unique_ptr<TH1> FDFlux_osc(
        dynamic_cast<TH1 *>(linsolv.FDFlux_osc->Clone()));
    FDFlux_osc->SetDirectory(nullptr);

    RelativeFluxMatches.emplace_back();
    for (int i = 0; i < Match->GetXaxis()->GetNbins(); ++i) {

      double rel_diff =
          (Match->GetBinContent(i + 1) - FDFlux_unosc->GetBinContent(i + 1)) /
          FDFlux_osc->GetBinContent(i + 1);

      RelativeFluxMatches.back().push_back(rel_diff);

      RelativeFluxDifferences->Fill(Match->GetXaxis()->GetBinCenter(i + 1),
                                    rel_diff);
    }

    if (SaveThrows) {
      Match->SetName((std::string("Match_") + std::to_string(t_it)).c_str());
      FDFlux_unosc->SetName(
          (std::string("Flux_Unosc_") + std::to_string(t_it)).c_str());
      FDFlux_osc->SetName(
          (std::string("Flux_osc") + std::to_string(t_it)).c_str());

      Match->SetDirectory(oupFile);
      FDFlux_unosc->SetDirectory(oupFile);
      FDFlux_osc->SetDirectory(oupFile);

      Match.release();
      FDFlux_unosc.release();
      FDFlux_osc.release();
    }
  }

  ThrowCoefficients->SetDirectory(oupFile);
  ThrowCoefficients.release();

  CovarianceBuilder cb(RelativeFluxMatches);

  std::unique_ptr<TMatrixD> covmat = GetTMatrixD(cb.GetCovMatrix());
  covmat->Write("difference_covmat");

  FillHistFromEigenVector(FDFlux_nom.get(), cb.GetMeanVector(), 0,
                          cb.GetStdDevVector());
  FDFlux_nom->SetName("difference_mean");
  FDFlux_nom->SetDirectory(oupFile);
  FDFlux_nom.release();

  RelativeFluxDifferences->SetDirectory(oupFile);
  RelativeFluxDifferences.release();

  oupFile->Write();
  oupFile->Close();
}
