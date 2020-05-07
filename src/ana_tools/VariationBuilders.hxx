#ifndef VARIATIONBUILDER_HXX_SEEN
#define VARIATIONBUILDER_HXX_SEEN

#include "PolyResponse.hxx"
#include "ROOTUtility.hxx"

#include "TDirectory.h"
#include "TH1.h"

#include "fhiclcpp/ParameterSet.h"

#include "Eigen/Dense"

#include <memory>
#include <string>
#include <vector>

class VariationBuilder {
protected:
  size_t NDOF;

  fhicl::ParameterSet paramset;

  TDirectory *diagdir;
  bool DumpDiagnostics;
  bool DumpMatrices;
  std::vector<std::vector<double>>
      diags_Predictions; // only used for diagnostics

  std::vector<std::string> Configurations;
  std::vector<std::string> Species;

  std::string Name;
  std::vector<double> NominalPrediction;
  std::vector<double> NominalPredictionError;

  std::vector<std::vector<double>> RelativeTweaks;
  std::vector<std::unique_ptr<TH1>> NominalHistogramSet;
  Eigen::MatrixXd CovarianceComponent;

  void BaseConfigure(fhicl::ParameterSet const &);

public:
  VariationBuilder() : diagdir(nullptr), DumpDiagnostics(false) {}
  virtual void Configure(fhicl::ParameterSet const &) = 0;
  virtual void Process() = 0;
  virtual std::vector<std::unique_ptr<TH1>> GetOneSigmaVariation() = 0;

  void SetDiagnosticDirectory(TDirectory *td) { diagdir = td; }
  std::vector<std::unique_ptr<TH1>>
  GetNominalHistograms(std::string const &suff = "") const {
    return CloneHistVector(NominalHistogramSet, suff);
  }
  Eigen::MatrixXd GetCovarianceComponent() const { return CovarianceComponent; }

  size_t GetNDOF() { return NDOF; }

  virtual ~VariationBuilder() = default;
};

class ThrownVariations : public VariationBuilder {

public:
  std::vector<std::unique_ptr<TH1>> GetOneSigmaVariation() {
    std::cout << "[ERROR]: Thrown variation (" << Name << ") cannot process directly to "
                 "OneSigmaVariations, Please use PCA tools instead."
              << std::endl;
    abort();
  }
  ThrownVariations() : VariationBuilder() {}
  void Configure(fhicl::ParameterSet const &);
  void Process();
};

class DiscreteVariations : public VariationBuilder {
  std::vector<double> sig_vals;
  std::vector<PolyResponse<2>> InterpolatedResponses;

  std::vector<std::vector<double>> diags_Errors; // only used for diagnostics

  bool fEvalAtOne;

public:
  DiscreteVariations() : VariationBuilder() {}
  std::vector<std::unique_ptr<TH1>> GetOneSigmaVariation();
  void Configure(fhicl::ParameterSet const &);
  void Process();
};

class DirectVariations : public VariationBuilder {
public:
  DirectVariations() : VariationBuilder() {}
  std::vector<std::unique_ptr<TH1>> GetOneSigmaVariation();
  void Configure(fhicl::ParameterSet const &);
  void Process();
};

class UniformVariations : public VariationBuilder {
public:
  UniformVariations() : VariationBuilder() {}
  std::vector<std::unique_ptr<TH1>> GetOneSigmaVariation();
  void Configure(fhicl::ParameterSet const &);
  void Process();
};

inline std::unique_ptr<VariationBuilder>
GetVariationBuilder(fhicl::ParameterSet const &ps) {
  std::unique_ptr<VariationBuilder> var;
  std::string const &type = ps.get<std::string>("Type", "");
  if (type == "Thrown") {
    var = std::unique_ptr<VariationBuilder>(new ThrownVariations());
  } else if (type == "Discrete") {
    var = std::unique_ptr<VariationBuilder>(new DiscreteVariations());
  } else if (type == "Direct") {
    var = std::unique_ptr<VariationBuilder>(new DirectVariations());
  } else if (type == "Uniform") {
    var = std::unique_ptr<VariationBuilder>(new UniformVariations());
  } else {
    std::cout << "[ERROR]: Invalid uncertainty component type, expected either "
                 "\"Thrown\" or \"Discrete\", but found: \""
              << type << "\"" << std::endl;
    throw;
  }
  var->Configure(ps);
  return var;
}

#endif
