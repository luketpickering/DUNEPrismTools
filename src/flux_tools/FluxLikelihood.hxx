#ifndef FLUXLIKELIHOOD_SEEN
#define FLUXLIKELIHOOD_SEEN

#include "TF1.h"
#include "TH1D.h"
#include "TDirectory.h"

#ifdef USE_FHICL_CPP
#include "fhiclcpp/ParameterSet.h"
#endif

#include <memory>
#include <string>
#include <vector>
#include <sstream>

class FluxLinearCombiner {
  std::vector<bool> ApplyRegularization;

protected:
  std::vector<std::unique_ptr<TH1D>> FluxSlices;
  mutable std::unique_ptr<TH1D> SummedFlux;

public:
  struct FluxSliceOptions {
    std::string InputFluxFile, InputFluxName;
    std::vector<std::pair<double, double>> XRanges;

    int MergeOAPBins;
    int MergeENuBins;

    double RegFactor;
  };

  FluxSliceOptions fFSO;

  void Initialize(FluxSliceOptions const &);
  void LoadFluxes();

  size_t GetNearestPeakingFluxIndex(double enu) const;
  virtual size_t GetNearestPeakingFluxIndex() const = 0;

  void SumFluxes(double const *FluxCoefficients) const;
  size_t GetNFluxes() const { return FluxSlices.size(); }
  double GetRegularizationPenalty(double const *FluxCoefficients) const;

  virtual double GetLikelihood(double const *coefficients) const = 0;

  virtual double GetPeakNorm() const = 0;

  virtual std::string State() const { return ""; };

  virtual void Write(TDirectory *) {};
};

class FluxTargetLikelihood : public FluxLinearCombiner {
  std::unique_ptr<TH1D> TargetFlux;
  std::unique_ptr<TH1D> InputFlux;

  Int_t fFitBinLow;
  Int_t fFitBinHigh;
  double fTargetPeakNorm;
  double fNDOverFDFitScaleFactor;

public:
  struct FluxTargetOptions {
    FluxSliceOptions input_flux_opts;

    /// If using a FitBetween mode:
    /// 0: Ignore all bins outside range
    /// 1: Try to force bins to 0
    /// 2: Exponential decay from target flux at closest kept bin.
    enum OutOfRangeModeEnum {
      kIgnore = 0,
      kZero,
      kExponentialDecay,
      kGaussianDecay
    };

    /// If using an OOR mode:
    /// 0: Include both out of ranges
    /// 1: Only include out of range to the left of the range
    /// 2: Only include out of range to the right of the range
    enum OutOfRangeSideEnum { kBoth = 0, kLeft, kRight };

    bool FitBetweenFoundPeaks;
    std::pair<double, double> FitBetween;

    int OutOfRangeMode;
    int OutOfRangeSide;
    double OORFactor;
    double ExpDecayRate;
    double TargetFractionalFreedom;

    bool UseNuPrismChi2;
  };

  struct LHood {
    double InFitRegionChi2;
    double OutOfFitRegionLH;
    double RegLH;
  };

  void Initialize(FluxTargetLikelihood::FluxTargetOptions const &);

  void SetTargetFlux(TH1D const *);
  using FluxLinearCombiner::GetNearestPeakingFluxIndex;
  size_t GetNearestPeakingFluxIndex() const;

  double BinDiff(Int_t bin_it) const;
  LHood GetLikelihoodComponents(double const *coefficients) const;

  mutable LHood last_lh;
  double GetLikelihood(double const *coefficients) const {
    last_lh = GetLikelihoodComponents(coefficients);
    return last_lh.InFitRegionChi2 + last_lh.OutOfFitRegionLH + last_lh.RegLH;
  }

  std::string State() const {
    std::stringstream ss("");
    ss << "LHood: { "
       << " InFitRegionChi2: " << last_lh.InFitRegionChi2
       << ", OutOfFitRegionLH: " << last_lh.OutOfFitRegionLH
       << ", RegLH: " << last_lh.RegLH;
    return ss.str();
  };

  double GetPeakNorm() const { return fTargetPeakNorm; }

  void Write(TDirectory *);
private:
  FluxTargetOptions fLHoodOpts;
};

class GausTargetLikelihood : public FluxLinearCombiner {
  std::unique_ptr<TF1> TargetGauss;

  Int_t fFitBinLow;
  Int_t fFitBinHigh;
  double fTargetPeakNorm;

public:
  struct GausTargetOptions {
    FluxSliceOptions input_flux_opts;

    double GaussC;
    double GaussW;
    std::pair<double, double> FitBetween;

    bool UseNuPrismChi2;
  };

  void Initialize(GausTargetLikelihood::GausTargetOptions const &);

  void SetGaussParameters(double GaussC, double GaussW);
  using FluxLinearCombiner::GetNearestPeakingFluxIndex;
  size_t GetNearestPeakingFluxIndex() const;

  double GetLikelihood(double const *coefficients) const;

  double GetPeakNorm() const { return fTargetPeakNorm; }

private:
  GausTargetOptions fLHoodOpts;
};

#ifdef USE_FHICL_CPP
FluxLinearCombiner::FluxSliceOptions
MakeFluxSliceOptions(fhicl::ParameterSet const &);
FluxTargetLikelihood::FluxTargetOptions
MakeFluxTargetOptions(fhicl::ParameterSet const &);
GausTargetLikelihood::GausTargetOptions
MakeGausTargetOptions(fhicl::ParameterSet const &);
#endif
FluxLinearCombiner::FluxSliceOptions MakeFluxSliceOptions(int argc,
                                                          char const *argv[]);
FluxTargetLikelihood::FluxTargetOptions
MakeFluxTargetOptions(int argc, char const *argv[]);
GausTargetLikelihood::GausTargetOptions
MakeGausTargetOptions(int argc, char const *argv[]);

#endif
