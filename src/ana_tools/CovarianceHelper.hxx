#ifndef COVARIANCEHELPER_HXX_SEEN
#define COVARIANCEHELPER_HXX_SEEN

#include "Eigen/Dense"

#include <memory>

class CovarianceBuilder {
  Eigen::MatrixXd CovMatrix;
  Eigen::VectorXd MeanVector;
  size_t NThrows_MeanCalc, NThrows_CovMatCalc;

  int NRows;
  bool MeanCalcFinalize;

  bool ZeroMean, MeanIsSet;

public:
  CovarianceBuilder(int NRows = 0);

  void Resize(int);

  void AddThrow_MeanCalc(double const *t);
  void FinalizeMeanCalc();

  void SetZeroMean();
  void SetMean(double const *t);

  void AddThrow_CovMatCalc(double const *t);
  void FinalizeCovMatCalc();

  Eigen::MatrixXd const &GetCovMatrix() { return CovMatrix; };
  Eigen::MatrixXd GetCorrMatrix();
};

#endif
