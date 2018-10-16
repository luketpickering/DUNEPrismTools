#include "CovarianceHelper.hxx"
#include "VALORModelClassifier.hxx"

#include <iomanip>
#include <limits>
#include <string>
#include <vector>

std::string oupfile = "TestMatrixThrows.root";

std::string inpfile;
std::string intmd;

Long64_t NThrows = 1;
Long64_t MatrixSize = std::numeric_limits<Long64_t>::max();
Long64_t NEvals = 20;

bool TestThrows = false;
bool AddOffDiag = false;

void SayUsage(char const *argv[]) {
  std::cout << "[USAGE]: " << argv[0]
            << "\n"
               "\t-o <friendoutput.root>     : TFile for output friend tree. \n"
               "\t-i <input.root,matrix>      :...\n"
               "\t-N <NThrows>     : Number of throws to make.\n"
               "\t--size           : Limit input matrix size.\n"
               "\t-E <NEigenvals>  : Ignore any more than -E evals."
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
    } else if (std::string(argv[opt]) == "-N") {
      NThrows = str2T<Long64_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-E") {
      NEvals = str2T<Long64_t>(argv[++opt]);
    } else if (std::string(argv[opt]) == "-i") {
      std::vector<std::string> params =
          ParseToVect<std::string>(argv[++opt], ",");
      if (params.size() != 2) {
        std::cout << "[ERROR]: Recieved " << params.size()
                  << " entrys for -i, expected 2." << std::endl;
        exit(1);
      }
      inpfile = params[0];
      intmd = params[1];
    } else if (std::string(argv[opt]) == "--size") {
      MatrixSize = str2T<Long64_t>(argv[++opt]);
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

  std::unique_ptr<TMatrixD> mat = GetROOTObject_uptr<TMatrixD>(inpfile, intmd);

  Eigen::MatrixXd UncertMatrix = GetEigenMatrix(mat.get());

  std::vector<std::vector<double>> ParameterThrows;

  ParameterThrows.resize(NThrows);

  if (UncertMatrix.rows() > MatrixSize) {
    Eigen::MatrixXd swap = UncertMatrix.block(0, 0, MatrixSize, MatrixSize);
    UncertMatrix = swap;
  } else if (UncertMatrix.rows() < MatrixSize) {
    MatrixSize = UncertMatrix.rows();
  }

  for (int t_it = 0; t_it < NThrows; ++t_it) {
    ParameterThrows[t_it].resize(MatrixSize);
  }

  EigenvalueHelper eh;
  eh.ComputeFromMatrix(UncertMatrix, true, NEvals);

  Eigen::MatrixXd ev_rebuild = Eigen::MatrixXd::Zero(MatrixSize, MatrixSize);

  Eigen::MatrixXd evectors = eh.GetEffectiveParameterVectors();

  for (int v_it = 0; v_it < evectors.cols(); ++v_it) {
    ev_rebuild += evectors.col(v_it) * evectors.col(v_it).transpose();
  }

  CovarianceThrower *ct_ev =
      new CovarianceEVThrower(UncertMatrix, true, NEvals);

  for (int t_it = 0; t_it < NThrows; ++t_it) {
    Eigen::VectorXd tvect = ct_ev->Throw();
    // #pragma omp parallel for
    for (int p_it = 0; p_it < MatrixSize; ++p_it) {
      ParameterThrows[t_it][p_it] = tvect(p_it);
    }

    if (t_it && !(t_it % 1000)) {
      std::cout << "[INFO]: Done " << t_it << "/" << NThrows << " throws."
                << std::endl;
    }
  }

  CovarianceBuilder cb_ev(ParameterThrows);

  TMatrixD DiffCovmat_EV(MatrixSize, MatrixSize);

  for (int i = 0; i < MatrixSize; i++) {
    for (int j = 0; j < MatrixSize; j++) {
      DiffCovmat_EV[i][j] = cb_ev.GetCovMatrix()(i, j) - UncertMatrix(i, j);
    }
  }

  TFile *of = new TFile(oupfile.c_str(), "RECREATE");

  TDirectory *tdir = of->mkdir("ThrowerTests");
  tdir->cd();

  std::unique_ptr<TMatrixD> input_covmat = GetTMatrixD(UncertMatrix);
  input_covmat->Write("input_covmat");

  std::unique_ptr<TMatrixD> rebuilt_covmat_EV =
      GetTMatrixD(cb_ev.GetCovMatrix());
  rebuilt_covmat_EV->Write("rebuilt_covmat_EV");

  std::unique_ptr<TMatrixD> rebuilt_covmat_EV_nothrow = GetTMatrixD(ev_rebuild);
  rebuilt_covmat_EV_nothrow->Write("rebuilt_covmat_EV_nothrow");

  DiffCovmat_EV.Clone()->Write("DiffCovmat_EV");

  of->Write();
  of->Close();
}
