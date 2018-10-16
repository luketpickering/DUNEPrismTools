#include "CovarianceHelper.hxx"

#include "Eigen/Cholesky"
#include "Eigen/Eigenvalues"

#include <SymEigsSolver.h>

#include <chrono>
#include <iostream>

Eigen::MatrixXd CovToCorr(Eigen::MatrixXd const &CovMatrix) {
  Eigen::MatrixXd CorrMat =
      Eigen::MatrixXd::Zero(CovMatrix.rows(), CovMatrix.rows());

  size_t const nr = CovMatrix.rows();
  for (size_t ix = 0; ix < nr; ++ix) {
    for (size_t jy = 0; jy < nr; ++jy) {

      CorrMat(ix, jy) = CovMatrix(ix, jy) /
                        (sqrt(CovMatrix(ix, ix)) * sqrt(CovMatrix(jy, jy)));
    }
  }
  return CorrMat;
}

CovarianceBuilder::CovarianceBuilder(int NRows)
    : CovMatrix(Eigen::MatrixXd::Zero(NRows, NRows)),
      MeanVector(Eigen::VectorXd::Zero(NRows)), NThrows_MeanCalc(0),
      NThrows_CovMatCalc(0), MeanCalcFinalize(false), ZeroMean(false),
      MeanIsSet(false), NRows(NRows) {}

CovarianceBuilder::CovarianceBuilder(
    std::vector<std::vector<double>> const &RandomVectors)
    : CovarianceBuilder(RandomVectors.size() ? RandomVectors.back().size()
                                             : 0) {
  if (!RandomVectors.size()) {
    return;
  }

  for (std::vector<double> const &flux_tweak : RandomVectors) {
    AddThrow_MeanCalc(flux_tweak.data());
  }

  FinalizeMeanCalc();

  for (std::vector<double> const &flux_tweak : RandomVectors) {
    AddThrow_CovMatCalc(flux_tweak.data());
  }
  FinalizeCovMatCalc();
}

void CovarianceBuilder::SetZeroMean() {
  ZeroMean = true;
  MeanVector = Eigen::VectorXd::Zero(NRows);
  MeanCalcFinalize = true;
}

void CovarianceBuilder::Resize(int NewNRows) {
  CovMatrix = Eigen::MatrixXd::Zero(NewNRows, NewNRows);
  MeanVector = Eigen::VectorXd::Zero(NewNRows);
  NThrows_MeanCalc = 0;
  NThrows_CovMatCalc = 0;
  NRows = NewNRows;
  MeanCalcFinalize = false;
  ZeroMean = false;
  MeanIsSet = false;
}

void CovarianceBuilder::AddThrow_MeanCalc(double const *t) {
  if (MeanCalcFinalize) {
    return;
  }

  NThrows_MeanCalc++;

  Eigen::Map<Eigen::VectorXd const> t_v(t, NRows);
  MeanVector += t_v;
}

void CovarianceBuilder::FinalizeMeanCalc() {
  if (MeanCalcFinalize || MeanIsSet) {
    return;
  }

  MeanVector.array() /= double(NThrows_MeanCalc);

  MeanCalcFinalize = true;
}

void CovarianceBuilder::SetMean(double const *t) {
  Eigen::Map<Eigen::VectorXd const> t_v(t, NRows);
  MeanVector = t_v;

  MeanIsSet = true;
}

void CovarianceBuilder::AddThrow_CovMatCalc(double const *t) {
  FinalizeMeanCalc();
  NThrows_CovMatCalc++;

  Eigen::Map<Eigen::VectorXd const> t_v(t, NRows);

  CovMatrix += (t_v - MeanVector) * (t_v - MeanVector).transpose();
}

void CovarianceBuilder::FinalizeCovMatCalc() {
  if ((!ZeroMean) && (!MeanIsSet) && (NThrows_MeanCalc != NThrows_CovMatCalc)) {
    std::cout << "[ERROR]: Added " << NThrows_MeanCalc
              << " throws to the mean calculator and " << NThrows_CovMatCalc
              << " to the covariance matrix builder." << std::endl;
    throw;
  }

  CovMatrix.bottomLeftCorner(NRows - 1, NRows - 1) =
      CovMatrix.topRightCorner(NRows - 1, NRows - 1).transpose();

  CovMatrix.array() /= (double(NThrows_CovMatCalc) - 1);
}

Eigen::MatrixXd CovarianceBuilder::GetCorrMatrix() {
  return CovToCorr(GetCovMatrix());
}

Eigen::VectorXd CovarianceBuilder::GetStdDevVector() {
  return CovMatrix.diagonal().cwiseSqrt();
}

void EigenvalueHelper::ComputeFromMatrix(Eigen::MatrixXd const &matrix,
                                         bool use_Spectra, size_t NEvals) {

  if (use_Spectra) {
    NEvals = std::min(NEvals, size_t(matrix.rows() - 1));

    Spectra::DenseSymMatProd<double> op(matrix);
    Spectra::SymEigsSolver<double, Spectra::LARGEST_ALGE,
                           Spectra::DenseSymMatProd<double>>
        eigs(&op, NEvals, std::min(2 * NEvals, size_t(matrix.rows())));
    // Initialize and compute
    eigs.init();
    eigs.compute();
    if (eigs.info() != Spectra::SUCCESSFUL) {
      std::cout << "[WARN]: Spectra Failed to find the top " << NEvals
                << " eigenvalues and vectors." << std::endl;
      return;
    }

    // Retrieve results
    EigVals = eigs.eigenvalues();
    EigVects = eigs.eigenvectors(NEvals);

    std::cout << "[INFO]: Spectra decomposition: EVect(" << EigVects.rows()
              << " x " << EigVects.cols() << ")" << std::endl;

  } else { // Use Eigen

    NEvals = std::min(NEvals, size_t(matrix.rows()));

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigs(matrix);

    // Retrieve results
    EigVals = eigs.eigenvalues();
    EigVects = eigs.eigenvectors();

    std::cout << "[INFO]: Eigen decomposition: EVect(" << EigVects.rows()
              << " x " << EigVects.cols() << ")" << std::endl;
  }

  std::vector<std::pair<double, size_t>> ev_ind;
  Eigen::MatrixXd EigVects_copy = EigVects;

  for (int i = 0; i < EigVals.rows(); ++i) {
    ev_ind.push_back({EigVals[i], i});
  }
  std::sort(
      ev_ind.begin(), ev_ind.end(),
      [](std::pair<double, size_t> const &l,
         std::pair<double, size_t> const &r) { return l.first > r.first; });

  for (int i = 0; i < EigVals.rows(); ++i) {
    EigVals[i] = ev_ind[i].first;
    EigVects.col(i) = EigVects_copy.col(ev_ind[i].second);
  }

  double sum_ev = 0;
  for (int i = 0; i < EigVals.rows(); ++i) {
    double mag2 = EigVects.col(i).dot(EigVects.col(i));
    std::cout << "[EigenDecomp]: EV[" << i << "]: " << EigVals(i)
              << ", Mag of evect = " << sqrt(mag2) << std::endl;
    sum_ev += EigVals(i);
  }
  std::cout << "[EigenDecomp]: Sum of kept eigen values = " << sum_ev
            << std::endl;
}

Eigen::MatrixXd EigenvalueHelper::GetEffectiveParameterVectors() {
  if (!EigVals.rows()) {
    return Eigen::MatrixXd();
  }
  Eigen::MatrixXd TweakParams = EigVects;

  for (int i = 0; i < EigVals.rows(); ++i) {
    TweakParams.col(i) *= sqrt(EigVals(i));
  }

  return TweakParams;
}

void CovarianceLDecompThrower::SetCovmat(Eigen::MatrixXd covmat) {
  Eigen::LLT<Eigen::MatrixXd> decomp(covmat.rows());
  decomp.compute(covmat);
  if (decomp.info() == Eigen::NumericalIssue) {
    throw std::runtime_error("Eigen failed to Cholesky decompose matrix.");
  }
  LMat = decomp.matrixL();
}

Eigen::VectorXd CovarianceLDecompThrower::Throw() {
  Eigen::VectorXd Throw(LMat.rows());
  for (int i = 0; i < LMat.rows(); ++i) {
    Throw(i) = (*RNJesus)(*RNEngine);
  }

  return LMat * Throw;
}

void CovarianceEVThrower::SetCovmat(Eigen::MatrixXd covmat) {
  eh.ComputeFromMatrix(covmat, UseSpectra, NEvals);
  ScaledEVects = eh.GetEffectiveParameterVectors();
}

Eigen::VectorXd CovarianceEVThrower::Throw() {
  int NVects = ScaledEVects.cols();

  Eigen::VectorXd rnd_vect(NVects);

  for (int p = 0; p < NVects; ++p) {
    rnd_vect(p) = (*RNJesus)(*RNEngine);
  }

  return ScaledEVects * rnd_vect;
}
