#include "CovarianceHelper.hxx"

#include "TDecompChol.h"

#include <iostream>

//#define COVHELPER_DEBUG

CovarianceBuilder::CovarianceBuilder(int NRows)
    : CovMatrix(Eigen::MatrixXd::Zero(NRows, NRows)),
      MeanVector(Eigen::VectorXd::Zero(NRows)), NThrows_MeanCalc(0),
      NThrows_CovMatCalc(0), NRows(NRows), MeanCalcFinalize(false),
      ZeroMean(false), MeanIsSet(false) {}

void CovarianceBuilder::SetZeroMean() {
  ZeroMean = true;
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

  #pragma omp parallel for
  for (int i = 0; i < NRows; ++i) {
#ifdef COVHELPER_DEBUG
    if (t[i] != t[i] || (t[i] && !std::isnormal(t[i]))) {
      std::cout << "[ERROR]: Attempted to add mean vector with bad value: "
                << t[i] << " at index " << i << std::endl;
      throw;
    }
#endif
    MeanVector(i) += t[i];
  }
}

void CovarianceBuilder::FinalizeMeanCalc() {
  if (MeanCalcFinalize || MeanIsSet) {
    return;
  }

#ifdef COVHELPER_DEBUG
  for (int i = 0; i < NRows; ++i) {

    if (MeanVector(i) != MeanVector(i) ||
        (MeanVector(i) && !std::isnormal(MeanVector(i)))) {
      std::cout << "[ERROR]: During mean finalizing, found bad value: "
                << MeanVector(i) << " at index " << i << std::endl;
      throw;
    }
  }
#endif

  double NThrowsWeight = 1.0 / double(NThrows_MeanCalc);

#ifdef COVHELPER_DEBUG
  if (!NThrowsWeight || !std::isnormal(NThrowsWeight)) {
    std::cout << "[ERROR]: Bad mean normalization weight: " << NThrowsWeight
              << std::endl;
    throw;
  }
#endif

  #pragma omp parallel for
  for (int i = 0; i < NRows; ++i) {
    MeanVector(i) *= NThrowsWeight;
  }

  MeanCalcFinalize = true;
}

void CovarianceBuilder::SetMean(double const *t) {

  for (int i = 0; i < NRows; ++i) {
#ifdef COVHELPER_DEBUG
    if (t[i] != t[i] || (t[i] && !std::isnormal(t[i]))) {
      std::cout << "[ERROR]: Attempted to set mean vector with bad value: "
                << t[i] << " at index " << i << std::endl;
      throw;
    }
#endif
    MeanVector(i) = t[i];
  }

  MeanIsSet = true;
}

void CovarianceBuilder::AddThrow_CovMatCalc(double const *t) {
  FinalizeMeanCalc();
  NThrows_CovMatCalc++;

#ifdef COVHELPER_DEBUG
  for (int i = 0; i < NRows; ++i) {

    if (t[i] != t[i] || (t[i] && !std::isnormal(t[i]))) {
      std::cout << "[ERROR]: Attempted to add throw vector with bad value: "
                << t[i] << " at index " << i << std::endl;
      throw;
    }
    if (!ZeroMean) {
      if (MeanVector(i) != MeanVector(i) ||
          (MeanVector(i) && !std::isnormal(MeanVector(i)))) {
        std::cout << "[ERROR]: Using non-zero mean and found bad mean value: "
                  << MeanVector(i) << " at index " << i << std::endl;
        throw;
      }
    }
  }
#endif

#ifdef COVHELPER_DEBUG
  for (int i = 0; i < NRows; ++i) {
    for (int j = 0; j < NRows; ++j) {

      if (CovMatrix(i, j) != CovMatrix(i, j) ||
          (CovMatrix(i, j) && !std::isnormal(CovMatrix(i, j)))) {
        std::cout << "[ERROR]: Found bad value in covariance matrix before "
                     "adding new throw ("
                  << NThrows_CovMatCalc << "): " << CovMatrix(i, j)
                  << " at index " << i << ", " << j << std::endl;
        throw;
      }
    }
  }
#endif

  if (ZeroMean) {
    #pragma omp parallel for
    for (int i = 0; i < NRows; ++i) {
      for (int j = 0; j < NRows; ++j) {
        CovMatrix(i, j) += t[i] * t[j];
      }
    }
  } else {

    #pragma omp parallel for
    for (int i = 0; i < NRows; ++i) {
      for (int j = 0; j < NRows; ++j) {
        CovMatrix(i, j) += (t[i] - MeanVector(i)) * (t[j] - MeanVector(j));
      }
    }
  }

#ifdef COVHELPER_DEBUG
  for (int i = 0; i < NRows; ++i) {
    for (int j = 0; j < NRows; ++j) {

      if (CovMatrix(i, j) != CovMatrix(i, j) ||
          (CovMatrix(i, j) && !std::isnormal(CovMatrix(i, j)))) {
        std::cout << "[ERROR]: Found bad value in covariance matrix "
                     "after new throw ("
                  << NThrows_CovMatCalc << "): " << CovMatrix(i, j)
                  << " at index " << i << ", " << j << std::endl;
        throw;
      }
    }
  }
#endif
}

void CovarianceBuilder::FinalizeCovMatCalc() {
  if ((!ZeroMean) && (!MeanIsSet) && (NThrows_MeanCalc != NThrows_CovMatCalc)) {
    std::cout << "[ERROR]: Added " << NThrows_MeanCalc
              << " throws to the mean calculator and " << NThrows_CovMatCalc
              << " to the covariance matrix builder." << std::endl;
    throw;
  }

#ifdef COVHELPER_DEBUG
  for (int i = 0; i < NRows; ++i) {
    for (int j = 0; j < NRows; ++j) {

      if (CovMatrix(i, j) != CovMatrix(i, j) ||
          (CovMatrix(i, j) && !std::isnormal(CovMatrix(i, j)))) {
        std::cout << "[ERROR]: Found bad value in covariance matrix: "
                  << CovMatrix(i, j) << " at index " << i << ", " << j
                  << std::endl;
        throw;
      }
    }
  }
#endif

  double NThrowsWeight = 1.0 / double(NThrows_CovMatCalc);

#ifdef COVHELPER_DEBUG
  if (!NThrowsWeight || !std::isnormal(NThrowsWeight)) {
    std::cout << "[ERROR]: Bad covariance throw normalization weight: "
              << NThrowsWeight << std::endl;
    throw;
  }
#endif

  #pragma omp parallel for
  for (int i = 0; i < NRows; ++i) {
    for (int j = 0; j < NRows; ++j) {
      CovMatrix(i, j) *= NThrowsWeight;
    }
  }
}

Eigen::MatrixXd CovarianceBuilder::GetCorrMatrix() {
  Eigen::MatrixXd CorrMat = Eigen::MatrixXd::Zero(NRows, NRows);

  #pragma omp parallel for
  for (Int_t ix = 0; ix < NRows; ++ix) {
    for (Int_t jy = 0; jy < NRows; ++jy) {

      CorrMat(ix, jy) = CovMatrix(ix, jy) /
                        (sqrt(CovMatrix(ix, ix)) * sqrt(CovMatrix(jy, jy)));
    }
  }
  return CorrMat;
}
