#ifndef EVRATELINEARSOLVER_STANDALONE_HXX_SEEN
#define EVRATELINEARSOLVER_STANDALONE_HXX_SEEN

#include "Eigen/Dense"

#include <iostream>

enum EvRateSolver { kSVD = 1, kQR, kNormal, kInverse };

Eigen::VectorXd SolveEvRate(Eigen::MatrixXd const &NDEvRateMatrix,
                            Eigen::VectorXd const &FDEvRate,
                            EvRateSolver algo_id = kInverse,
                            double regfactor = 1E-10) {

  bool use_reg = regfactor > 0;
  size_t NFluxes = NDEvRateMatrix.cols();
  size_t NBins = NDEvRateMatrix.rows();
  size_t NEqs = NBins + (use_reg * NFluxes);

  Eigen::MatrixXd NDEvRateMatrix_wreg;
  Eigen::VectorXd FDEvRate_wreg;

  if (use_reg) {
    NDEvRateMatrix_wreg = Eigen::MatrixXd::Zero(NEqs, NFluxes);
    FDEvRate_wreg = Eigen::VectorXd::Zero(NEqs);
    FDEvRate_wreg.topRows(NBins) = FDEvRate.topRows(NBins);

    for (size_t row_it = 0; row_it < (NFluxes - 1); ++row_it) {
      NDEvRateMatrix_wreg(row_it + NBins, row_it) = regfactor;
    }
    NDEvRateMatrix_wreg(NEqs - 1, NFluxes - 1) = regfactor;
  }

  Eigen::VectorXd solution;

  switch (algo_id) {
  case kSVD: {
    if (use_reg) {
      solution =
          NDEvRateMatrix_wreg.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
              .solve(FDEvRate_wreg);
    } else {
      solution =
          NDEvRateMatrix.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
              .solve(FDEvRate);
    }
    break;
  }
  case kQR: {
    if (use_reg) {
      solution = NDEvRateMatrix_wreg.colPivHouseholderQr().solve(FDEvRate_wreg);
    } else {
      solution = NDEvRateMatrix.colPivHouseholderQr().solve(FDEvRate);
    }
    break;
  }
  case kNormal: {
    if (use_reg) {
      solution = (NDEvRateMatrix_wreg.transpose() * NDEvRateMatrix_wreg)
                     .ldlt()
                     .solve(NDEvRateMatrix_wreg.transpose() * FDEvRate_wreg);
    } else {
      solution = (NDEvRateMatrix.transpose() * NDEvRateMatrix)
                     .topRows(NBins)
                     .ldlt()
                     .solve(NDEvRateMatrix.transpose() * FDEvRate);
    }
    break;
  }
  case kInverse: {
    if (use_reg) {
      solution = ((NDEvRateMatrix_wreg.topRows(NBins).transpose() *
                   NDEvRateMatrix_wreg.topRows(NBins)) +
                  NDEvRateMatrix_wreg.bottomRows(NFluxes).transpose() *
                      NDEvRateMatrix_wreg.bottomRows(NFluxes))
                     .inverse() *
                 NDEvRateMatrix_wreg.topRows(NBins).transpose() *
                 FDEvRate_wreg.topRows(NBins);
    } else {
      solution = (NDEvRateMatrix.transpose() * NDEvRateMatrix).inverse() *
                 NDEvRateMatrix.transpose() * FDEvRate;
    }
    break;
  }
  }

  if (!solution.rows()) {
    std::cout << "[ERROR]: Failed to perform Event rate linear algebra."
              << std::endl;
    exit(0);
  }
  return solution;
}
#endif
