#include "FluxFitter.hxx"

#include "SliceConfigTreeReader.hxx"

#include "ROOTUtility.hxx"
#include "StringParserUtility.hxx"

#include "Fit/Fitter.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

#include <algorithm>
#include <chrono>
#include <functional>

//#define DEBUG_LHOOD

void FluxFitter::Initialize(FluxFitterOptions const &ffo) { fFitterOpts = ffo; }

void FluxFitter::InitializeFlux(
    FluxFitter::FluxFitterOptions const &ffo,
    FluxTargetLikelihood::FluxTargetOptions const &fto) {
  Initialize(ffo);
  fFitterOpts.IsGauss = false;
  std::unique_ptr<FluxTargetLikelihood> fluxLHood =
      std::make_unique<FluxTargetLikelihood>();

  fluxLHood->Initialize(fto);
  LHoodEval = std::move(fluxLHood);

  if (!fFitterOpts.InputTargetFluxFile.size() ||
      !fFitterOpts.InputTargetFluxName.size()) {
    std::cout << "[ERROR]: No target flux file specified." << std::endl;
    throw;
  }

  tgtFlux = GetHistogram_uptr<TH1D>(fFitterOpts.InputTargetFluxFile,
                                    fFitterOpts.InputTargetFluxName);

  SetTargetFlux(tgtFlux.get(), true);
}
void FluxFitter::InitializeGauss(
    FluxFitter::FluxFitterOptions const &ffo,
    GausTargetLikelihood::GausTargetOptions const &gto) {
  Initialize(ffo);
  fFitterOpts.IsGauss = true;
  std::unique_ptr<GausTargetLikelihood> gausLHood =
      std::make_unique<GausTargetLikelihood>();

  gausLHood->Initialize(gto);

  LHoodEval = std::move(gausLHood);

  InitializeCoefficientGuess();
}

void FluxFitter::InitializeCoefficientGuess() {

  Coefficients.clear();
  size_t NFluxes = LHoodEval->GetNFluxes();
  std::fill_n(std::back_inserter(Coefficients), NFluxes, 0);

  if (fFitterOpts.InpCoeffFile.size()) {
    if (fFitterOpts.InpCoeffDir.size() &&
        (fFitterOpts.InpCoeffDir.back() != '/')) {
      fFitterOpts.InpCoeffDir += "/";
    }

    std::cout << "[INFO]: Reading input previous fit results." << std::endl;
    SliceConfig sc(fFitterOpts.InpCoeffFile, fFitterOpts.InpCoeffDir);

    std::vector<std::pair<double, double>> inpXRanges = sc.GetXRanges();
    std::vector<double> inpCoeffs = sc.GetCoeffs();

    if (inpXRanges.size() != LHoodEval->fFSO.XRanges.size()) {
      std::cout << "[ERROR]: Input fit result had " << inpXRanges.size()
                << " off-axis slices, but here we have found "
                << LHoodEval->fFSO.XRanges.size()
                << ". Input results are incompatible." << std::endl;
      throw;
    }

    for (size_t f_it = 0; f_it < inpXRanges.size(); ++f_it) {
      if ((fabs(inpXRanges[f_it].first -
                (100. * LHoodEval->fFSO.XRanges[f_it].first)) > 1E-5) ||
          (fabs(inpXRanges[f_it].second -
                (100. * LHoodEval->fFSO.XRanges[f_it].second)) > 1E-5)) {
        std::cout << "\tFlux window[" << f_it << "] = {"
                  << inpXRanges[f_it].first << " -- " << inpXRanges[f_it].second
                  << "}, Coeff = " << inpCoeffs[f_it] << std::endl;
        std::cout << "[ERROR]: Here we found Flux window[" << f_it << "] = {"
                  << LHoodEval->fFSO.XRanges[f_it].first << " -- "
                  << LHoodEval->fFSO.XRanges[f_it].second
                  << "}. Input results are incompatible." << std::endl;
        throw;
      }
      Coefficients[f_it] = inpCoeffs[f_it];
    }

  } else {
    size_t flux_with_closest_peak_it = LHoodEval->GetNearestPeakingFluxIndex();

    for (size_t flux_it = 0; flux_it < NFluxes; flux_it++) {
      if (abs(flux_it - flux_with_closest_peak_it) ==
          1) { // Neighbours are negative
        Coefficients[flux_it] = -0.35;
      } else if (abs(flux_it - flux_with_closest_peak_it) == 0) {
        Coefficients[flux_it] = 1;
      } else { // Others start free
        Coefficients[flux_it] = 0;
      }
    }
  }
}

bool FluxFitter::Fit() {

  size_t NFluxes = LHoodEval->GetNFluxes();

  if (!fMinimizer) {
    fMinimizer = std::unique_ptr<ROOT::Math::Minimizer>(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  }

  fMinimizer->SetMaxFunctionCalls(fFitterOpts.MaxLikelihoodCalls);
  fMinimizer->SetTolerance(fFitterOpts.MINUITTolerance);

  for (size_t flux_it = 0; flux_it < NFluxes; flux_it++) {
    std::stringstream ss("");
    ss << "fluxpar" << flux_it;

    fMinimizer->SetVariable(flux_it, ss.str().c_str(), Coefficients[flux_it],
                            0.1);
    fMinimizer->SetVariableLimits(flux_it, -fFitterOpts.CoeffLimit,
                                  fFitterOpts.CoeffLimit);
  }

  size_t fcn_call = 0;
  auto lap = std::chrono::high_resolution_clock::now();
  ROOT::Math::Functor fcn(
      [&](double const *coeffs) {
        if (fcn_call && !(fcn_call % 5000)) {
          auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(
              std::chrono::high_resolution_clock::now() - lap);
          lap = std::chrono::high_resolution_clock::now();
          std::cout << "[INFO]: Call: " << fcn_call << ", Last 5k calls took "
                    << (double(diff.count()) / 1000.0) << "s" << std::endl;
          std::cout << "[INFO]: Likelihood state: " << LHoodEval->State()
                    << std::endl;
        }
#ifdef DEBUG_LHOOD
        std::cout << "[INFO]: Step " << fcn_call << std::endl;
        for (size_t c = 0; c < NFluxes; ++c) {
          std::cout << "\t c[" << c << "] = " << coeffs[c] << std::endl;
        }
#endif
        double lh = LHoodEval->GetLikelihood(coeffs);
#ifdef DEBUG_LHOOD
        std::cout << "LH = " << lh << std::endl;
#endif

        fcn_call++;
        return lh;
      },
      NFluxes);

  fMinimizer->SetFunction(fcn);

  int fitstatus = fMinimizer->Minimize();
  if (fitstatus) {
    std::cout << "[WARN]: Failed to find minimum (STATUS: " << fitstatus << ")."
              << std::endl;
  }
  fMinimizer->PrintResults();

  for (size_t flux_it = 0; flux_it < NFluxes; flux_it++) {
    Coefficients[flux_it] = fMinimizer->X()[flux_it];
  }

  return (!fitstatus);
}

void FluxFitter::SetCoefficients(std::vector<double> coeffs_new) {
  Coefficients = std::move(coeffs_new);
}
void FluxFitter::SetTargetFlux(TH1D const *tgt, bool GuessCoefficients) {
  if (fFitterOpts.IsGauss) {
    std::cout << "[ERROR]: Attempted to SetTargetFlux on a FluxFitter that is "
                 "configured to fit gaussian fluxes."
              << std::endl;
    throw;
  }

  dynamic_cast<FluxTargetLikelihood *>(LHoodEval.get())->SetTargetFlux(tgt);

  if (GuessCoefficients) {
    InitializeCoefficientGuess();
  }
}
void FluxFitter::SetGaussParameters(double GaussC, double GaussW,
                                    bool GuessCoefficients) {
  if (!fFitterOpts.IsGauss) {
    std::cout << "[ERROR]: Attempted to SetGaussParameters on a FluxFitter "
                 "that is configured to fit target fluxes."
              << std::endl;
    throw;
  }
  dynamic_cast<GausTargetLikelihood *>(LHoodEval.get())
      ->SetGaussParameters(GaussC, GaussW);

  if (GuessCoefficients) {
    InitializeCoefficientGuess();
  }
}

#ifdef USE_FHICL_CPP
FluxFitter::FluxFitterOptions
MakeFluxFitterOptions(fhicl::ParameterSet const &);
#endif
FluxFitter::FluxFitterOptions MakeFluxFitterOptions(int argc,
                                                    char const *argv[]) {

  FluxFitter::FluxFitterOptions ffo;

  ffo.IsGauss = false;
  ffo.InpCoeffFile = "";
  ffo.InpCoeffDir = "";
  ffo.InputTargetFluxFile = "";
  ffo.InputTargetFluxName = "";
  ffo.CoeffLimit = 10;
  ffo.MaxLikelihoodCalls = 5E5;
  ffo.MINUITTolerance = 1E-5;

  int opt = 1;
  while (opt < argc) {
    if (std::string(argv[opt]) == "-A") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      ffo.InpCoeffFile = params[0];
      if (params.size() > 1) {
        ffo.InpCoeffDir = params[1];
      } else {
        ffo.InpCoeffDir = "";
      }
    } else if (std::string(argv[opt]) == "-t") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -t, expected 2." << std::endl;
        throw;
      }
      ffo.InputTargetFluxFile = params[0];
      ffo.InputTargetFluxName = params[1];
    } else if (std::string(argv[opt]) == "-n") {
      ffo.MaxLikelihoodCalls = str2T<int>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-T") {
      ffo.MINUITTolerance = str2T<double>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-c") {
      ffo.CoeffLimit = str2T<double>(argv[++opt]);
    }
    opt++;
  }
  return ffo;
}
