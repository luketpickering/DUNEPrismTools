#ifndef COVARIANCEHELPER_HXX_SEEN
#define COVARIANCEHELPER_HXX_SEEN

#include "Eigen/Dense"

#include "TMatrixD.h"

#include "ROOTUtility.hxx"

#include <memory>
#include <random>
#include <vector>

Eigen::MatrixXd CovToCorr(Eigen::MatrixXd const &);

class CovarianceBuilder {
  Eigen::MatrixXd CovMatrix;
  Eigen::VectorXd MeanVector;
  size_t NThrows_MeanCalc, NThrows_CovMatCalc;

  bool MeanCalcFinalize;

  bool ZeroMean, MeanIsSet;

public:
  int NRows;

  CovarianceBuilder(int NRows = 0);
  CovarianceBuilder(std::vector<std::vector<double>> const &,
                    bool force_zero_mean = false);

  void Resize(int);

  void AddThrow_MeanCalc(double const *t);
  void FinalizeMeanCalc();

  void SetZeroMean();
  void SetMean(double const *t);

  void AddThrow_CovMatCalc(double const *t);
  void FinalizeCovMatCalc();

  Eigen::MatrixXd const &GetCovMatrix() { return CovMatrix; };
  Eigen::VectorXd const &GetMeanVector() { return MeanVector; };
  Eigen::VectorXd GetStdDevVector();
  Eigen::MatrixXd GetCorrMatrix();
};

class EigenvalueHelper {
public:
  Eigen::VectorXd EigVals;
  Eigen::MatrixXd EigVects;

  void ComputeFromMatrix(Eigen::MatrixXd const &, bool use_Spectra = true,
                         size_t NEvals = 50);

  Eigen::MatrixXd GetEffectiveParameterVectors();
};

class CovarianceThrower {
protected:
  std::unique_ptr<std::mt19937> RNEngine;
  std::unique_ptr<std::normal_distribution<double>> RNJesus;

public:
  CovarianceThrower() {
    std::random_device r;
    RNEngine = std::unique_ptr<std::mt19937>(new std::mt19937(r()));
    RNJesus = std::unique_ptr<std::normal_distribution<double>>(
        new std::normal_distribution<double>(0, 1));
  }

  virtual void SetCovmat(Eigen::MatrixXd) = 0;
  virtual Eigen::VectorXd Throw() = 0;
};

class CovarianceLDecompThrower : public CovarianceThrower {
  Eigen::MatrixXd LMat;

public:
  CovarianceLDecompThrower(Eigen::MatrixXd covmat) : CovarianceThrower() {
    SetCovmat(std::move(covmat));
  }

  void SetCovmat(Eigen::MatrixXd);
  Eigen::VectorXd Throw();
};

class CovarianceEVThrower : public CovarianceThrower {
  EigenvalueHelper eh;
  Eigen::MatrixXd ScaledEVects;

  bool UseSpectra;
  size_t NEvals;

public:
  CovarianceEVThrower(Eigen::MatrixXd covmat, bool us, size_t ne)
      : CovarianceThrower(), UseSpectra(us), NEvals(ne) {
    SetCovmat(std::move(covmat));
  }

  void SetCovmat(Eigen::MatrixXd);
  Eigen::VectorXd Throw();
};

#endif
