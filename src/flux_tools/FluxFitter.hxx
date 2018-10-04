#ifndef FLUXFITTER_SEEN
#define FLUXFITTER_SEEN

#include "FluxLikelihood.hxx"

#include "TDirectory.h"

#include "Math/Minimizer.h"

#include <memory>
#include <string>
#include <vector>

class FluxFitter {
  std::unique_ptr<FluxLinearCombiner> LHoodEval;
  std::unique_ptr<TH1D> tgtFlux;

  std::vector<double> Coefficients;
  std::unique_ptr<ROOT::Math::Minimizer> fMinimizer;

public:
  struct FluxFitterOptions {
    bool IsGauss;

    std::string InpCoeffFile, InpCoeffDir;

    std::string InputTargetFluxFile, InputTargetFluxName;

    double CoeffLimit;

    double MINUITTolerance;
    unsigned int MaxLikelihoodCalls;
  };

  void InitializeFlux(FluxFitterOptions const &,
                      FluxTargetLikelihood::FluxTargetOptions const &);
  void InitializeGauss(FluxFitterOptions const &,
                       GausTargetLikelihood::GausTargetOptions const &);

  bool Fit();

  void SetCoefficients(std::vector<double>);
  void SetTargetFlux(TH1D const *, bool GuessCoefficients = false);
  void SetGaussParameters(double GaussC, double GaussW,
                          bool GuessCoefficients = false);

  void Write(TDirectory *td){
    LHoodEval->Write(td);
  }

private:
  void Initialize(FluxFitterOptions const &);
  void InitializeCoefficientGuess();

  FluxFitterOptions fFitterOpts;
};

#ifdef USE_FHICL_CPP
FluxFitter::FluxFitterOptions
MakeFluxFitterOptions(fhicl::ParameterSet const &);
#endif
FluxFitter::FluxFitterOptions MakeFluxFitterOptions(int argc,
                                                    char const *argv[]);

#endif
