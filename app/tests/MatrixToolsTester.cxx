#include "CovarianceHelper.h"
#include "EDepTreeReader.h"
#include "VALORModelClassifier.h"

#include "Utils.hxx"

#include "omp.h"

#include <iomanip>
#include <string>
#include <vector>

std::string oupfile = "TestMatrixThrows.root";

std::string CovmatFile, CovmatHist;
Long64_t NThrows = 1;

bool TestThrows = false;
bool AddOffDiag = false;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-o <friendoutput.root>     : TFile for output friend tree. \n"
               "\t-C <covmatfile.root,hname> : Input covariance matrix to "
               "throw xsec weights from.\n"
               "\t-N <NThrows>     : Number of throws to make.\n"
            << std::endl;
}

void handleOpts(int argc, char const *argv[]) {
  int opt = 1;
  while (opt < argc) {
    if ((std::string(argv[opt]) == "-?") ||
        (std::string(argv[opt]) == "--help")) {
      SayUsage(argv);
      exit(0);
    } else if (std::string(argv[opt]) == "-o") {
      oupfile = argv[++opt];
    } else if (std::string(argv[opt]) == "-C") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -C, expected 2." << std::endl;
        throw;
      }
      CovmatFile = params[0];
      CovmatHist = params[1];
    } else if (std::string(argv[opt]) == "-N") {
      NThrows = str2T<Long64_t>(argv[++opt]);
    } else {
      std::cout << "[ERROR]: Unknown option: " << argv[opt] << std::endl;
      SayUsage(argv);
      exit(1);
    }
    opt++;
  }
}

int main(int argc, char const *argv[]) {
  TH1::SetDefaultSumw2();
  handleOpts(argc, argv);

  TMatrixDSym *UncertMatrix;
  if (CovmatFile.size() && CovmatHist.size()) {
  } else {
    UncertMatrix =
        new TMatrixDSym(static_cast<int>(VALORModel::TrueClass::kNVALORDials));
    for (int d_it = 0;
         d_it < static_cast<int>(VALORModel::TrueClass::kNVALORDials); ++d_it) {
      (*UncertMatrix)[d_it][d_it] =
          pow(GetDefaultDialTweak(static_cast<VALORModel::TrueClass>(d_it)), 2);
    }
    if (AddOffDiag) {
      for (int d_it = 0;
           d_it < (static_cast<int>(VALORModel::TrueClass::kNVALORDials) - 1);
           ++d_it) {
        (*UncertMatrix)[d_it][d_it + 1] = (*UncertMatrix)[d_it][d_it] * 0.1;
        (*UncertMatrix)[d_it + 1][d_it] = (*UncertMatrix)[d_it][d_it] * 0.1;
      }
    }
  }
  std::vector<std::vector<double> > ParameterThrows;

  ParameterThrows.resize(NThrows);
  for (int t_it = 0; t_it < NThrows; ++t_it) {
    ParameterThrows[t_it].resize(
        static_cast<int>(VALORModel::TrueClass::kNVALORDials));
  }
  CovarianceThrower ct(*UncertMatrix);

  for (int t_it = 0; t_it < NThrows; ++t_it) {
    TMatrixD const *CovVector = ct.Throw();
    for (int p_it = 0;
         p_it < static_cast<int>(VALORModel::TrueClass::kNVALORDials); ++p_it) {
      ParameterThrows[t_it][p_it] = (*CovVector)[p_it][0];
    }
  }

  CovarianceBuilder cb(static_cast<int>(VALORModel::TrueClass::kNVALORDials));

  for (int t_it = 0; t_it < NThrows; ++t_it) {
    cb.AddThrow_MeanCalc(ParameterThrows[t_it].data());
  }

  for (int t_it = 0; t_it < NThrows; ++t_it) {
    cb.AddThrow_CovMatCalc(ParameterThrows[t_it].data());
  }

  cb.FinalizeCovMatCalc();

  TMatrixD DiffCovmat(static_cast<int>(VALORModel::TrueClass::kNVALORDials),
                      static_cast<int>(VALORModel::TrueClass::kNVALORDials));

  for (int i = 0; i < static_cast<int>(VALORModel::TrueClass::kNVALORDials);
       i++) {
    for (int j = 0; j < static_cast<int>(VALORModel::TrueClass::kNVALORDials);
         j++) {
      DiffCovmat[i][j] = ((*cb.GetCovMatrix())[i][j] - (*UncertMatrix)[i][j]);
    }
  }

  TFile *of = new TFile(oupfile.c_str(), "RECREATE");

  TDirectory *tdir = of->mkdir("ThrowerTests");
  tdir->cd();

  ct.Write();
  cb.Write();
  DiffCovmat.Clone()->Write("DifferenceCovMats");

  of->cd();
  TDirectory *edir = of->mkdir("EigenValueTests");
  edir->cd();

  TMatrixD Var1(10, 1);
  TMatrixD Var2(10, 1);
  TMatrixD Var3(10, 1);

  Var1[0][0] = 1;
  Var1[1][0] = 1;
  Var1[2][0] = 1;
  Var1[3][0] = 1;
  Var1[4][0] = 1;
  Var1[5][0] = -1;
  Var1[6][0] = -1;
  Var1[7][0] = -1;
  Var1[8][0] = -1;
  Var1[9][0] = -1;

  Var2[0][0] = 1;
  Var2[1][0] = 1;
  Var2[2][0] = -1;
  Var2[3][0] = -1;
  Var2[4][0] = -1;
  Var2[5][0] = -1;
  Var2[6][0] = -1;
  Var2[7][0] = -1;
  Var2[8][0] = 1;
  Var2[9][0] = 1;

  Var3[0][0] = 1;
  Var3[1][0] = 2;
  Var3[2][0] = 3;
  Var3[3][0] = 4;
  Var3[4][0] = 5;
  Var3[5][0] = 5;
  Var3[6][0] = 4;
  Var3[7][0] = 3;
  Var3[8][0] = 2;
  Var3[9][0] = 1;

  CovarianceBuilder cb2(10);
  cb2.SetZeroMean();

  std::cout << "[INFO]: Input variations:" << std::endl;

  std::cout.precision(3);

  std::cout << "\tV1: {" << std::flush;
  for (int i = 0; i < 10; ++i) {
    std::cout << Var1[i][0] << ((i == 9) ? " }" : ", ") << std::flush;
  }
  std::cout << std::endl;
  std::cout << "\tV2: {" << std::flush;
  for (int i = 0; i < 10; ++i) {
    std::cout << Var2[i][0] << ((i == 9) ? " }" : ", ") << std::flush;
  }
  std::cout << std::endl;
  std::cout << "\tV3: {" << std::flush;
  for (int i = 0; i < 10; ++i) {
    std::cout << Var3[i][0] << ((i == 9) ? " }" : ", ") << std::flush;
  }
  std::cout << std::endl;

  TRandom3 rnjesus;
  TMatrixD Resp(10, 1);
  for (int t_it = 0; t_it < NThrows; ++t_it) {
    double p1 = rnjesus.Gaus();
    double p2 = rnjesus.Gaus();
    double p3 = rnjesus.Gaus();

    for (int i = 0; i < 10; ++i) {
      Var1[i][0] *= p1;
      Var2[i][0] *= p2;
      Var3[i][0] *= p3;
    }

    Resp = Var1 + Var2 + Var3;

    cb2.AddThrow_CovMatCalc(&Resp);

    for (int i = 0; i < 10; ++i) {
      Var1[i][0] /= p1;
      Var2[i][0] /= p2;
      Var3[i][0] /= p3;
    }
  }

  cb2.FinalizeCovMatCalc();

  TAxis ax(10, 0, 10);
  EVCovMatWeightEngine ev(cb2.GetCovMatrix(), &ax);

  TMatrixD evmat(*ev.GetEVMatrix());
  TVectorD evvect(*ev.GetEVVect());

  cb2.Write();
  ev.Write();

  of->cd();
  TDirectory *etdir = of->mkdir("EigenValueReThrow");
  etdir->cd();

  for (int i = 0; i < 10; ++i) {
    Var1[i][0] = evmat[i][0];
    Var2[i][0] = evmat[i][1];
    Var3[i][0] = evmat[i][2];
  }

  std::cout << "[INFO]: EV variations:" << std::endl;

  std::cout.precision(3);

  std::cout << "\tV1: {" << std::flush;
  for (int i = 0; i < 10; ++i) {
    std::cout << Var1[i][0] << ((i == 9) ? " }" : ", ") << std::flush;
  }
  std::cout << std::endl;
  std::cout << "\tV2: {" << std::flush;
  for (int i = 0; i < 10; ++i) {
    std::cout << Var2[i][0] << ((i == 9) ? " }" : ", ") << std::flush;
  }
  std::cout << std::endl;
  std::cout << "\tV3: {" << std::flush;
  for (int i = 0; i < 10; ++i) {
    std::cout << Var3[i][0] << ((i == 9) ? " }" : ", ") << std::flush;
  }
  std::cout << std::endl;

  CovarianceBuilder cb3(10);
  cb3.SetZeroMean();

  // double sumev = evvect[0] + evvect[1] + evvect[2];
  for (int t_it = 0; t_it < NThrows; ++t_it) {
    double p1 =
        rnjesus.Gaus() * sqrt(fabs(evvect[0])) * (evvect[0] > 0 ? 1 : -1);
    double p2 =
        rnjesus.Gaus() * sqrt(fabs(evvect[1])) * (evvect[1] > 0 ? 1 : -1);
    double p3 =
        rnjesus.Gaus() * sqrt(fabs(evvect[2])) * (evvect[2] > 0 ? 1 : -1);

    for (int i = 0; i < 10; ++i) {
      Var1[i][0] *= p1;
      Var2[i][0] *= p2;
      Var3[i][0] *= p3;
    }

    Resp = Var1 + Var2 + Var3;

    cb3.AddThrow_CovMatCalc(&Resp);

    for (int i = 0; i < 10; ++i) {
      Var1[i][0] /= p1;
      Var2[i][0] /= p2;
      Var3[i][0] /= p3;
    }
  }

  cb3.FinalizeCovMatCalc();

  EVCovMatWeightEngine ev2(cb3.GetCovMatrix(), &ax);

  TMatrixD ReThrowDiff(10, 10);

  for (int i = 0; i < 10; i++) {
    for (int j = 0; j < 10; j++) {
      ReThrowDiff[i][j] =
          ((*cb3.GetCovMatrix())[i][j] - (*cb2.GetCovMatrix())[i][j]);
    }
  }
  ReThrowDiff.Clone()->Write("ReThrowDifferenceCovMats");

  cb3.Write();

  of->Write();
  of->Close();
}
