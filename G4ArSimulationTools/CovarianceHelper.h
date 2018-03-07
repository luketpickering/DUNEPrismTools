#include "TAxis.h"
#include "TDecompChol.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "TVectorD.h"

#include <iostream>

class CovarianceBuilder {
  TMatrixD *CovMatrix;
  TMatrixD *MeanVector;
  size_t NThrows_MeanCalc, NThrows_CovMatCalc;

  int NRows;
  bool MeanCalcFinalize;

  bool ZeroMean;

 public:
  CovarianceBuilder(int NRows)
      : NThrows_MeanCalc(0),
        NThrows_CovMatCalc(0),
        NRows(NRows),
        MeanCalcFinalize(false),
        ZeroMean(false) {
    CovMatrix = new TMatrixD(NRows, NRows);
    MeanVector = new TMatrixD(NRows, 1);

    for (int i = 0; i < NRows; ++i) {
      (*MeanVector)[i][0] = 0;
      for (int j = 0; j < NRows; ++j) {
        (*CovMatrix)[i][j] = 0;
      }
    }
  }

  void SetZeroMean() {
    ZeroMean = true;
    MeanCalcFinalize = true;
  }

  void AddThrow_MeanCalc(double *t) {
    TMatrixD v(NRows, 1);
    for (int i = 0; i < NRows; ++i) {
      v[i][0] = t[i];
    }

    AddThrow_MeanCalc(&v);
  }

  void AddThrow_MeanCalc(TMatrixD const *Vector) {
    if (MeanCalcFinalize) {
      return;
    }

    NThrows_MeanCalc++;
    for (int i = 0; i < NRows; ++i) {
      (*MeanVector)[i][0] += (*Vector)[i][0];
    }
  }

  void FinalizeMeanCalc() {
    if (MeanCalcFinalize) {
      return;
    }

    double NThrowsWeight = 1.0 / double(NThrows_MeanCalc);

    for (int i = 0; i < NRows; ++i) {
      (*MeanVector)[i][0] *= NThrowsWeight;
    }

    MeanCalcFinalize = true;
  }

  void AddThrow_CovMatCalc(double *t) {
    TMatrixD v(NRows, 1);
    for (int i = 0; i < NRows; ++i) {
      v[i][0] = t[i];
    }

    AddThrow_CovMatCalc(&v);
  }

  void AddThrow_CovMatCalc(TMatrixD const *Vector) {
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

  void FinalizeCovMatCalc() {
    if (!ZeroMean && (NThrows_MeanCalc != NThrows_CovMatCalc)) {
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

  TMatrixD *GetCovMatrix() { return CovMatrix; }

  void Write() {
    CovMatrix->Write("CovMatVariations");
    if (!ZeroMean) {
      MeanVector->Write("MeanVariations");
    }
  }
};

class CovarianceThrower {
  TMatrixD *UncertMatrix;
  TMatrixD *LMatrix;
  TMatrixD *RVector;

  TMatrixD *CVector;

  TRandom3 *RNJesus;

  int NRows;

  CovarianceThrower(int NRows) : NRows(NRows) {
    UncertMatrix = new TMatrixD(NRows, NRows);

    RNJesus = new TRandom3();

    LMatrix = new TMatrixD(NRows, NRows);
    RVector = new TMatrixD(NRows, 1);
    CVector = new TMatrixD(NRows, 1);
  }

 public:
  void SetupDecomp() {
    TDecompChol decomp(*UncertMatrix);
    if (!decomp.Decompose()) {
      std::cout << "[ERROR]: Failed to decompose uncertainty matrix."
                << std::endl;
      exit(1);
    }
    (*LMatrix) = decomp.GetU();
    (*LMatrix) = LMatrix->Transpose(*LMatrix);
  }

  CovarianceThrower(TMatrixDSym &covmat)
      : CovarianceThrower(covmat.GetNrows()) {
    for (int i = 0; i < NRows; ++i) {
      for (int j = 0; j < NRows; ++j) {
        (*UncertMatrix)[i][j] = covmat[i][j];
      }
    }

    SetupDecomp();
  }

  CovarianceThrower(TMatrixD &covmat) : CovarianceThrower(covmat.GetNrows()) {
    for (int i = 0; i < NRows; ++i) {
      for (int j = 0; j < NRows; ++j) {
        (*UncertMatrix)[i][j] = covmat[i][j];
      }
    }

    SetupDecomp();
  }

  TMatrixD const *Throw() {
    for (int p_it = 0; p_it < NRows; ++p_it) {
      (*RVector)[p_it][0] = RNJesus->Gaus();
    }
    (*CVector) = (*LMatrix) * (*RVector);

    return CVector;
  }

  void Write() {
    UncertMatrix->Write("UncertMatrix");
    LMatrix->Write("LMatrix");
  }
};

class EVCovMatWeightEngine {
  Int_t NRows;

  TMatrixD *UncertMatrix;
  TAxis *Axis;

  TMatrixD *EigenMatrix;
  TVectorD *EigenValues;

  void SetupDecomp() {
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

  EVCovMatWeightEngine(int NRows) : NRows(NRows) {
    UncertMatrix = new TMatrixD(NRows, NRows);
  }

 public:
  EVCovMatWeightEngine(TMatrixDSym const *UncertMatrix, TAxis const *Axis)
      : EVCovMatWeightEngine(UncertMatrix->GetNrows()) {
    for (int i = 0; i < NRows; ++i) {
      for (int j = 0; j < NRows; ++j) {
        (*this->UncertMatrix)[i][j] = (*UncertMatrix)[i][j];
      }
    }
    this->Axis = new TAxis(*Axis);

    SetupDecomp();
  }

  EVCovMatWeightEngine(TMatrixD const *UncertMatrix, TAxis const *Axis)
      : EVCovMatWeightEngine(UncertMatrix->GetNrows()) {
    for (int i = 0; i < NRows; ++i) {
      for (int j = 0; j < NRows; ++j) {
        (*this->UncertMatrix)[i][j] = (*UncertMatrix)[i][j];
      }
    }
    this->Axis = new TAxis(*Axis);

    SetupDecomp();
  }

  double GetWeight(int Bin_Id, int Param_id, double Param_Val) {
    return Param_Val * (*EigenValues)[Param_id] *
           (*EigenMatrix)[Bin_Id][Param_id];
  }
  double GetWeight(double BinnedVariable, int Param_id, double Param_Val) {
    return GetWeight(Axis->FindFixBin(BinnedVariable - 1), Param_id, Param_Val);
  }

  TMatrixD const *GetEVMatrix() { return EigenMatrix; }
  TVectorD const *GetEVVect() { return EigenValues; }

  void Write() {
    EigenMatrix->Write("DecompEigenVectors");
    EigenValues->Write("DecompEigenValues");

    for (int i = 0; i < NRows; ++i) {
      for (int j = 0; j < NRows; ++j) {
        (*this->EigenMatrix)[j][i] *= (*EigenValues)[i];
      }
    }

    EigenMatrix->Write("DecompEigenVectors_valnorm");
  }
};

void NormColumnMatrix(TMatrixD &UncertMatrix) {
  double sum = 0;
  for (Int_t i = 0; i < UncertMatrix.GetNrows(); ++i) {
    sum += (UncertMatrix[i][0]*UncertMatrix[i][0]);
  }
  for (Int_t i = 0; i < UncertMatrix.GetNrows(); ++i) {
    UncertMatrix[i][0] /= sqrt(sum);
  }
}
void NormColumnMatrix(TMatrixD *UncertMatrix) {
  NormColumnMatrix(*UncertMatrix);
}
