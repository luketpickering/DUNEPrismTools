#ifndef EVRATELINEARSOLVER_STANDALONE_HXX_SEEN
#define EVRATELINEARSOLVER_STANDALONE_HXX_SEEN

#include "Eigen/Dense"

enum EvRateSolver { kSVD = 1, kQR, kNormal, kInverse };

Eigen::VectorXd SolveEvRate(Eigen::MatrixXd const &EvRateMatrix,
                            Eigen::VectorXd const &Target,
                            EvRateSolver algo_id = kInverse,
                            double regfactor = 1E-10) {

  bool use_reg = reg_param > 0;
  size_t NFluxes = EvRateMatrix.cols();
  size_t NBins = EvRateMatrix.rows();
  size_t NEqs = NBins + (use_reg * NFluxes);

  Eigen::MatrixXd FluxMatrix_Solve = Eigen::MatrixXd::Zero(NEqs, NFluxes);
  Eigen::VectorXd Target_Solve = Eigen::MatrixXd::Zero(NEqs);
  Target_Solve.topRows(NBins) = Target.topRows(NBins);

  Eigen::VectorXd solution;

  if (use_reg) {
    for (size_t row_it = 0; row_it < (NFluxes - 1); ++row_it) {
      FluxMatrix_Solve(row_it + NBins, row_it) = reg_param;
    }
    FluxMatrix_Solve(NEqs - 1, NFluxes - 1) = reg_param;
  }

  switch (algo_id) {
  case kSVD: {
    if (use_reg) {
      solution =
          FluxMatrix_Solve.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
              .solve(Target_Solve);
    } else {
      solution = FluxMatrix_Solve.topRows(NBins)
                     .bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV)
                     .solve(Target_Solve.topRows(NBins));
    }
    break;
  }
  case kQR: {
    if (use_reg) {
      solution = FluxMatrix_Solve.colPivHouseholderQr().solve(Target_Solve);
    } else {
      solution = FluxMatrix_Solve.topRows(NBins).colPivHouseholderQr().solve(
          Target_Solve.topRows(NBins));
    }
    break;
  }
  case kNormal: {
    if (use_reg) {
      solution = (FluxMatrix_Solve.transpose() * FluxMatrix_Solve)
                     .ldlt()
                     .solve(FluxMatrix_Solve.transpose() * Target_Solve);
    } else {
      solution = (FluxMatrix_Solve.transpose() * FluxMatrix_Solve)
                     .topRows(NBins)
                     .ldlt()
                     .solve(FluxMatrix_Solve.topRows(NBins).transpose() *
                            Target_Solve.topRows(NBins));
    }
    break;
  }
  case kInverse: {
    if (use_reg) {
      solution = ((FluxMatrix_Solve.topRows(NBins).transpose() *
                   FluxMatrix_Solve.topRows(NBins)) +
                  FluxMatrix_Solve.bottomRows(NFluxes).transpose() *
                      FluxMatrix_Solve.bottomRows(NFluxes))
                     .inverse() *
                 FluxMatrix_Solve.topRows(NBins).transpose() *
                 Target_Solve.topRows(NBins);
    } else {
      solution = (FluxMatrix_Solve.topRows(NBins).transpose() *
                  FluxMatrix_Solve.topRows(NBins))
                     .inverse() *
                 FluxMatrix_Solve.topRows(NBins).transpose() *
                 Target_Solve.topRows(NBins);
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
