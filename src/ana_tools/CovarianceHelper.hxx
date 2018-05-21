#ifndef COVARIANCEHELPER_HXX_SEEN
#define COVARIANCEHELPER_HXX_SEEN

#include "TAxis.h"
#include "TH1.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "TVectorD.h"

class CovarianceBuilder {
  TMatrixD *CovMatrix;
  TMatrixD *CorrMatrix;
  TMatrixD *MeanVector;
  size_t NThrows_MeanCalc, NThrows_CovMatCalc;

  int NRows;
  bool MeanCalcFinalize;

  bool ZeroMean, MeanIsSet;

public:
  CovarianceBuilder(int NRows);

  void SetZeroMean();

  void AddThrow_MeanCalc(double *t);
  void AddThrow_MeanCalc(TH1 *t);
  void AddThrow_MeanCalc(TMatrixD const *Vector);

  void FinalizeMeanCalc();

  void SetMean(double *t);
  void SetMean(TH1 *t);
  void SetMean(TMatrixD const *Vector);

  void AddThrow_CovMatCalc(double *t);
  void AddThrow_CovMatCalc(TH1 *t);
  void AddThrow_CovMatCalc(TMatrixD const *Vector);

  void FinalizeCovMatCalc();

  TMatrixD *GetCovMatrix();
  TMatrixD *GetCorrMatrix();

  void Write();
};

class CovarianceThrower {
  TMatrixD *UncertMatrix;
  TMatrixD *LMatrix;
  TMatrixD *RVector;

  TMatrixD *CVector;

  TRandom3 *RNJesus;

  int NRows;

  UInt_t Seed;

  CovarianceThrower(int NRows, UInt_t Seed = 0);

public:
  void SetupDecomp(double decompTol = 0xdeadbeef);

  CovarianceThrower(TMatrixDSym &covmat, UInt_t Seed = 0,
                    double decompTol = 0xdeadbeef);
  CovarianceThrower(TMatrixD &covmat, UInt_t Seed = 0,
                    double decompTol = 0xdeadbeef);

  TMatrixD const *Throw();

  void Write();
};

class EVCovMatWeightEngine {
  Int_t NRows;

  TMatrixD *UncertMatrix;
  TAxis *Axis;

  TMatrixD *EigenMatrix;
  TVectorD *EigenValues;

  void SetupDecomp();

  EVCovMatWeightEngine(int NRows);

public:
  EVCovMatWeightEngine(TMatrixDSym const *UncertMatrix, TAxis const *Axis);

  EVCovMatWeightEngine(TMatrixD const *UncertMatrix, TAxis const *Axis);

  double GetWeight(int Bin_Id, int Param_id, double Param_Val);
  double GetWeight(double BinnedVariable, int Param_id, double Param_Val);

  TMatrixD const *GetEVMatrix();
  TVectorD const *GetEVVect();

  void Write();
};

void NormColumnMatrix(TMatrixD &UncertMatrix);
void NormColumnMatrix(TMatrixD *UncertMatrix);

#endif
