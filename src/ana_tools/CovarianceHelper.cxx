#include "CovarianceHelper.hxx"

#include "TDecompChol.h"

#include <iostream>

CovarianceBuilder::CovarianceBuilder(int NRows)
    : NThrows_MeanCalc(0),
      NThrows_CovMatCalc(0),
      NRows(NRows),
      MeanCalcFinalize(false),
      ZeroMean(false),
      MeanIsSet(false) {
  CovMatrix = new TMatrixD(NRows, NRows);
  CorrMatrix = nullptr;
  MeanVector = new TMatrixD(NRows, 1);

  for (int i = 0; i < NRows; ++i) {
    (*MeanVector)[i][0] = 0;
    for (int j = 0; j < NRows; ++j) {
      (*CovMatrix)[i][j] = 0;
    }
  }
}

void CovarianceBuilder::SetZeroMean() {
  ZeroMean = true;
  MeanCalcFinalize = true;
}

void CovarianceBuilder::AddThrow_MeanCalc(double *t) {
  TMatrixD v(NRows, 1);
  for (int i = 0; i < NRows; ++i) {
    v[i][0] = t[i];
  }

  AddThrow_MeanCalc(&v);
}

void CovarianceBuilder::AddThrow_MeanCalc(TH1 *t) {
  TMatrixD v(NRows, 1);
  for (int i = 0; i < NRows; ++i) {
    v[i][0] = t->GetBinContent(i + 1);
  }

  AddThrow_MeanCalc(&v);
}

void CovarianceBuilder::AddThrow_MeanCalc(TMatrixD const *Vector) {
  if (MeanCalcFinalize) {
    return;
  }

  NThrows_MeanCalc++;
  for (int i = 0; i < NRows; ++i) {
    (*MeanVector)[i][0] += (*Vector)[i][0];
  }
}

void CovarianceBuilder::FinalizeMeanCalc() {
  if (MeanCalcFinalize) {
    return;
  }

  double NThrowsWeight = 1.0 / double(NThrows_MeanCalc);

  for (int i = 0; i < NRows; ++i) {
    (*MeanVector)[i][0] *= NThrowsWeight;
  }

  MeanCalcFinalize = true;
}

void CovarianceBuilder::SetMean(double *t) {
  TMatrixD v(NRows, 1);
  for (int i = 0; i < NRows; ++i) {
    v[i][0] = t[i];
  }

  SetMean(&v);
}

void CovarianceBuilder::SetMean(TH1 *t) {
  TMatrixD v(NRows, 1);
  for (int i = 0; i < NRows; ++i) {
    v[i][0] = t->GetBinContent(i + 1);
  }

  SetMean(&v);
}

void CovarianceBuilder::SetMean(TMatrixD const *Vector) {
  for (int i = 0; i < NRows; ++i) {
    (*MeanVector)[i][0] = (*Vector)[i][0];
  }

  MeanIsSet = true;
}

void CovarianceBuilder::AddThrow_CovMatCalc(double *t) {
  TMatrixD v(NRows, 1);
  for (int i = 0; i < NRows; ++i) {
    v[i][0] = t[i];
  }

  AddThrow_CovMatCalc(&v);
}

void CovarianceBuilder::AddThrow_CovMatCalc(TH1 *t) {
  TMatrixD v(NRows, 1);
  for (int i = 0; i < NRows; ++i) {
    v[i][0] = t->GetBinContent(i + 1);
  }

  AddThrow_CovMatCalc(&v);
}

void CovarianceBuilder::AddThrow_CovMatCalc(TMatrixD const *Vector) {
  FinalizeMeanCalc();
  NThrows_CovMatCalc++;

  if (ZeroMean) {
    for (int i = 0; i < NRows; ++i) {
      for (int j = 0; j < NRows; ++j) {
        (*CovMatrix)[i][j] += (*Vector)[i][0] * (*Vector)[j][0];
      }
    }
  } else {
    for (int i = 0; i < NRows; ++i) {
      for (int j = 0; j < NRows; ++j) {
        (*CovMatrix)[i][j] += ((*Vector)[i][0] - (*MeanVector)[i][0]) *
                              ((*Vector)[j][0] - (*MeanVector)[j][0]);
      }
    }
  }
}

void CovarianceBuilder::FinalizeCovMatCalc() {
  if ((!ZeroMean) && (!MeanIsSet) &&
      (NThrows_MeanCalc != NThrows_CovMatCalc)) {
    std::cout << "[ERROR]: Added " << NThrows_MeanCalc
              << " throws to the mean calculator and " << NThrows_CovMatCalc
              << " to the covariance matrix builder." << std::endl;
    throw;
  }

  double NThrowsWeight = 1.0 / double(NThrows_CovMatCalc);

  for (int i = 0; i < NRows; ++i) {
    for (int j = 0; j < NRows; ++j) {
      (*CovMatrix)[i][j] *= NThrowsWeight;
    }
  }
}

TMatrixD *CovarianceBuilder::GetCovMatrix() { return CovMatrix; }

TMatrixD *CovarianceBuilder::GetCorrMatrix() {
  if (CorrMatrix) {
    delete CorrMatrix;
  }
  CorrMatrix = new TMatrixD(NRows, NRows);

  for (Int_t ix = 0; ix < NRows; ++ix) {
    for (Int_t jy = 0; jy < NRows; ++jy) {
      Double_t BinC = (*CovMatrix)[ix][jy];

      Double_t BinCx = (*CovMatrix)[ix][ix];
      Double_t BinCy = (*CovMatrix)[jy][jy];
      Double_t CVar = BinC / (sqrt(BinCx) * sqrt(BinCy));

      (*CorrMatrix)[ix][jy] = CVar;
    }
  }

  return CorrMatrix;
}

void CovarianceBuilder::Write() {
  CovMatrix->Write("CovMatVariations");
  if (!ZeroMean) {
    MeanVector->Write("MeanVariations");
  }
}



CovarianceThrower::CovarianceThrower(int NRows, UInt_t Seed) : NRows(NRows) {
  UncertMatrix = new TMatrixD(NRows, NRows);

  this->Seed = Seed;
  RNJesus = new TRandom3(Seed);

  LMatrix = new TMatrixD(NRows, NRows);
  RVector = new TMatrixD(NRows, 1);
  CVector = new TMatrixD(NRows, 1);
}

void CovarianceThrower::SetupDecomp() {
  TDecompChol decomp(*UncertMatrix);
  if (!decomp.Decompose()) {
    std::cout << "[ERROR]: Failed to decompose uncertainty matrix."
              << std::endl;
    exit(1);
  }
  (*LMatrix) = decomp.GetU();
  (*LMatrix) = LMatrix->Transpose(*LMatrix);
}

CovarianceThrower::CovarianceThrower(TMatrixDSym &covmat, UInt_t Seed)
    : CovarianceThrower(covmat.GetNrows(), Seed) {
  for (int i = 0; i < NRows; ++i) {
    for (int j = 0; j < NRows; ++j) {
      (*UncertMatrix)[i][j] = covmat[i][j];
    }
  }

  SetupDecomp();
}

CovarianceThrower::CovarianceThrower(TMatrixD &covmat, UInt_t Seed)
    : CovarianceThrower(covmat.GetNrows(), Seed) {
  for (int i = 0; i < NRows; ++i) {
    for (int j = 0; j < NRows; ++j) {
      (*UncertMatrix)[i][j] = covmat[i][j];
    }
  }

  SetupDecomp();
}

TMatrixD const *CovarianceThrower::Throw() {
  for (int p_it = 0; p_it < NRows; ++p_it) {
    (*RVector)[p_it][0] = RNJesus->Gaus();
  }
  (*CVector) = (*LMatrix) * (*RVector);

  return CVector;
}

void CovarianceThrower::Write() {
  UncertMatrix->Write("UncertMatrix");
  LMatrix->Write("LMatrix");
}

void EVCovMatWeightEngine::SetupDecomp() {
  EigenValues = new TVectorD(0);
  EigenMatrix = new TMatrixD(UncertMatrix->EigenVectors(*EigenValues));

  std::cout << "[INFO]: EV decomp: " << std::endl;

  for (Int_t i = 0; i < NRows; ++i) {
    std::cout << "\t[" << i << "] EVal: " << (*EigenValues)[i] << ", EVec: {"
              << std::flush;
    for (Int_t j = 0; j < NRows; ++j) {
      std::cout << (*EigenMatrix)[j][i] << ((j != (NRows - 1)) ? ", " : " }")
                << std::flush;
    }
    std::cout << std::endl;
  }
}

EVCovMatWeightEngine::EVCovMatWeightEngine(int NRows) : NRows(NRows) {
  UncertMatrix = new TMatrixD(NRows, NRows);
}

EVCovMatWeightEngine::EVCovMatWeightEngine(TMatrixDSym const *UncertMatrix, TAxis const *Axis)
    : EVCovMatWeightEngine(UncertMatrix->GetNrows()) {
  for (int i = 0; i < NRows; ++i) {
    for (int j = 0; j < NRows; ++j) {
      (*this->UncertMatrix)[i][j] = (*UncertMatrix)[i][j];
    }
  }
  this->Axis = new TAxis(*Axis);

  SetupDecomp();
}

EVCovMatWeightEngine::EVCovMatWeightEngine(TMatrixD const *UncertMatrix, TAxis const *Axis)
    : EVCovMatWeightEngine(UncertMatrix->GetNrows()) {
  for (int i = 0; i < NRows; ++i) {
    for (int j = 0; j < NRows; ++j) {
      (*this->UncertMatrix)[i][j] = (*UncertMatrix)[i][j];
    }
  }
  this->Axis = new TAxis(*Axis);

  SetupDecomp();
}

double EVCovMatWeightEngine::GetWeight(int Bin_Id, int Param_id, double Param_Val) {
  return Param_Val * (*EigenValues)[Param_id] *
         (*EigenMatrix)[Bin_Id][Param_id];
}
double EVCovMatWeightEngine::GetWeight(double BinnedVariable, int Param_id, double Param_Val) {
  return GetWeight(Axis->FindFixBin(BinnedVariable - 1), Param_id, Param_Val);
}

TMatrixD const *EVCovMatWeightEngine::GetEVMatrix() { return EigenMatrix; }
TVectorD const *EVCovMatWeightEngine::GetEVVect() { return EigenValues; }

void EVCovMatWeightEngine::Write() {
  EigenMatrix->Write("DecompEigenVectors");
  EigenValues->Write("DecompEigenValues");

  for (int i = 0; i < NRows; ++i) {
    for (int j = 0; j < NRows; ++j) {
      (*this->EigenMatrix)[j][i] *= (*EigenValues)[i];
    }
  }

  EigenMatrix->Write("DecompEigenVectors_valnorm");
}

void NormColumnMatrix(TMatrixD &UncertMatrix) {
  double sum = 0;
  for (Int_t i = 0; i < UncertMatrix.GetNrows(); ++i) {
    sum += (UncertMatrix[i][0] * UncertMatrix[i][0]);
  }
  for (Int_t i = 0; i < UncertMatrix.GetNrows(); ++i) {
    UncertMatrix[i][0] /= sqrt(sum);
  }
}
void NormColumnMatrix(TMatrixD *UncertMatrix) {
  NormColumnMatrix(*UncertMatrix);
}
